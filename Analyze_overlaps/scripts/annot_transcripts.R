#This script categorizes TT-seq transcripts as protein-coding, lncRNA, uasRNA and eRNA
library(tidyverse)
library(egg)
library(gridExtra)
library(viridis)
library(GenomicRanges)


theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), 
                  legend.text = element_text(size = 6),
                  axis.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.title.x = element_blank()))


day_vec <- c("d0", "d2", "d4")

args <- commandArgs(trailingOnly = TRUE)
files_dir <- args[1]
output_dir <- args[2]

annot_out <- lapply(day_vec, function(x)
{
gencode_plus <- read.delim(paste0(files_dir, "XX_", x, "_plus_transcripts_list.bed"), 
                           col.names = c("chr", "start", "end", "ensembl_id", "gencode_type", "id"), header = FALSE)   
gencode_minus <- read.delim(paste0(files_dir, "XX_", x, "_minus_transcripts_list.bed"), 
                           col.names = c("chr", "start", "end", "ensembl_id", "gencode_type", "id"), header = FALSE)


gencode_df <- rbind(gencode_plus, gencode_minus) %>% 
  select(id, ensembl_id, gencode_type)

assembly <- read.delim(paste0(files_dir, "XX_", x, "_trimmed_transcripts.bed"), 
                       col.names = c("chr", "start", "end", "id", "score", "strand"), header = FALSE)
source <-  read.delim(paste0(files_dir, x, "_XX_transcripts_source.txt"), 
                      header  = TRUE)

transcript_df <- left_join(assembly, source) %>% 
  left_join(gencode_df) %>% 
  na.omit()

#Get out lncRNA genes and write to GRanges
lnc_df <- transcript_df %>% 
  filter(!str_detect(gencode_type, ".*protein.*"))

lnc_gr <- GRanges(seqnames = lnc_df[[1]],
                  IRanges(start = lnc_df[[2]], end = lnc_df[[3]]), 
                  strand = lnc_df[[6]])

mcols(lnc_gr) <- DataFrame(id = lnc_df[[4]], chromhmm = lnc_df[[8]], ensembl_id = lnc_df[[10]],
                           gencode_type = lnc_df[[11]])

#Get out coding genes and write to gRanges
coding_df <- transcript_df %>% 
  filter(str_detect(gencode_type, ".*protein.*"))

coding_gr <- GRanges(seqnames = coding_df[[1]],
                  IRanges(start = coding_df[[2]], end = coding_df[[3]]), 
                  strand = coding_df[[6]])

mcols(coding_gr) <- DataFrame(id = coding_df[[4]], chromhmm = coding_df[[8]], ensembl_id = coding_df[[10]],
                           gencode_type = coding_df[[11]])

#Check if transcript is uaRNA (starting from opposite strand within 2 kb).
upstream_coding <- flank(coding_gr, width = 1000, start = TRUE)
strand(upstream_coding) <- ifelse(strand(upstream_coding) == "+", "-", "+")

ua_check <- findOverlaps(lnc_gr, upstream_coding, type = "start", maxgap = 500)
ua_gr <- lnc_gr[queryHits(ua_check)]

#Check if transcript is an enhancer RNA
eRNA_gr <- lnc_gr[mcols(lnc_gr)$chromhmm == "enhancer"]

#Check if lncRNAs overlap coding genes
genic_check <- findOverlaps(lnc_gr, coding_gr, ignore.strand = TRUE)
genic_gr <- lnc_gr[queryHits(genic_check)]
intergenic_gr <- lnc_gr[-queryHits(genic_check)]

#Assign types back to original dataframe
annotated_transcripts <- transcript_df %>% 
  mutate(annot = ifelse(id %in% coding_df$id, "protein-coding",
         ifelse(id %in% ua_gr$id, "uaRNA",
                ifelse(id %in% eRNA_gr$id, "eRNA",
                ifelse(id %in% genic_gr$id, "genic lncRNA", "intragenic lncRNA"))))) %>% 
  mutate(day = x)

write_delim(annotated_transcripts, paste0(output_dir, "XX_", x, "_full_annotation.txt"), delim = "\t")

return(annotated_transcripts)
})

out_df <- bind_rows(annot_out)

order <- rev(c("protein-coding", "genic lncRNA", "intragenic lncRNA", "eRNA", "uaRNA"))

plot_annot <- out_df %>% 
  ggplot() +
  geom_bar(aes(x = day, group = factor(annot, levels = order), fill = factor(annot, levels = order)), 
           color = NA, stat = "count", position = "stack") +
  geom_text(stat='count', aes(label=..count.., x = day, y = ..count.. + 4000), size = 6/2.8, color = "black") +
  labs(y = "Number of Transcripts") +
  ylim(0, NA) +
  scale_fill_viridis(discrete = TRUE, option = "C", begin = 0.1, end = 0.9)


pdf(paste0(output_dir, "FigS2E_TU_types_bar.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(plot_annot, width = unit(2.5, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))
dev.off()
  
