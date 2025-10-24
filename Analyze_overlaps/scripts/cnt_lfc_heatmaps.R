#Script to plot LFC heatmaps of Cut&Tag
library(bsseq)
library(tidyverse)
library(egg)
library(gridExtra)
library(viridis)
library(EnvStats)
library(Rsubread)
library(genomation)

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(),
                  axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
                  axis.text = element_text(size = 6), strip.text = element_text(size = 6)))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
bam_dir <- args[2]
fig_dir <- args[3]

day_list <- c("d0", "d2", "d4")

for (day in day_list ) 
{
#Reads in case and control files
case_3prime <- read.delim(paste0(data_dir, "XX_", day, "_3prime_case.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand"))
cont_3prime <- read.delim(paste0(data_dir, "XX_", day, "_3prime_control.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand"))

case_5prime <- read.delim(paste0(data_dir, "XX_", day, "_5prime_case.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand"))
cont_5prime <- read.delim(paste0(data_dir, "XX_", day, "_5prime_control.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand"))

#Reads in Info about corresponding cases and controls
case_controls_3prime <- read.delim(paste0(data_dir, "XX_", day, "_3prime_case_controls.txt")) %>% 
  select(case_id, control_id)
case_controls_5prime <- read.delim(paste0(data_dir, "XX_", day, "_5prime_case_controls.txt")) %>% 
  select(case_id, control_id)

#Limits promoter region to 500bp upstream and 200bp downstream
case_3prime_saf <- case_3prime %>%
  mutate(prom_start = ifelse(strand == "+", start - 500, end - 200),
         prom_end = ifelse(strand == "+", start + 200, end + 500)) %>% 
  transmute(GeneID = id, Chr = chr, Start = prom_start, End = prom_end, Strand = ".")

cont_3prime_saf <- cont_3prime %>%
  mutate(prom_start = ifelse(strand == "+", start - 500, end - 200),
         prom_end = ifelse(strand == "+", start + 200, end + 500)) %>% 
  transmute(GeneID = id, Chr = chr, Start = prom_start, End = prom_end, Strand = ".")

case_5prime_saf <- case_5prime %>%
  mutate(prom_start = ifelse(strand == "+", start - 500, end - 200),
         prom_end = ifelse(strand == "+", start + 200, end + 500)) %>% 
  transmute(GeneID = id, Chr = chr, Start = prom_start, End = prom_end, Strand = ".")

cont_5prime_saf <- cont_5prime %>%
  mutate(prom_start = ifelse(strand == "+", start - 500, end - 200),
         prom_end = ifelse(strand == "+", start + 200, end + 500)) %>% 
  transmute(GeneID = id, Chr = chr, Start = prom_start, End = prom_end, Strand = ".")

#Counts reads for cnt data in promoter range
setwd(bam_dir)
temp = normalizePath(list.files(pattern=paste0("XX_H.*", day, ".*.bam$"), full.names = TRUE))

case_3prime_featureCounts <- featureCounts(temp, annot.ext = case_3prime_saf, isPairedEnd = TRUE)
cont_3prime_featureCounts <- featureCounts(temp, annot.ext = cont_3prime_saf, isPairedEnd = TRUE)
case_5prime_featureCounts <- featureCounts(temp, annot.ext = case_5prime_saf, isPairedEnd = TRUE)
cont_5prime_featureCounts <- featureCounts(temp, annot.ext = cont_5prime_saf, isPairedEnd = TRUE)


#Calculates CPM for all promoters
total_reads_case_3prime <- case_3prime_featureCounts$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_case_3prime  <- sweep(as.matrix(case_3prime_featureCounts$counts), 2, total_reads_case_3prime, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000})  %>% 
  rownames_to_column("case_id") %>% 
  mutate(type = "3prime") %>% 
  left_join(case_controls_3prime)

names(counts_case_3prime) <- gsub(x = names(counts_case_3prime), pattern = "\\.bam", replacement = "")
names(counts_case_3prime) <- gsub(x = names(counts_case_3prime), pattern = "XX_", replacement = "case_")

total_reads_cont_3prime <- cont_3prime_featureCounts$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_cont_3prime  <- sweep(as.matrix(cont_3prime_featureCounts$counts), 2, total_reads_cont_3prime, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000})  %>% 
  rownames_to_column("control_id") %>% 
  mutate(type = "3prime") %>% 
  left_join(case_controls_3prime)

names(counts_cont_3prime) <- gsub(x = names(counts_cont_3prime), pattern = "\\.bam", replacement = "")
names(counts_cont_3prime) <- gsub(x = names(counts_cont_3prime), pattern = "XX_", replacement = "control_")

total_reads_case_5prime <- case_5prime_featureCounts$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_case_5prime  <- sweep(as.matrix(case_5prime_featureCounts$counts), 2, total_reads_case_5prime, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000})  %>% 
  rownames_to_column("case_id") %>% 
  mutate(type = "5prime") %>% 
  left_join(case_controls_5prime)

names(counts_case_5prime) <- gsub(x = names(counts_case_5prime), pattern = "\\.bam", replacement = "")
names(counts_case_5prime) <- gsub(x = names(counts_case_5prime), pattern = "XX_", replacement = "case_")

total_reads_cont_5prime <- cont_5prime_featureCounts$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_cont_5prime  <- sweep(as.matrix(cont_5prime_featureCounts$counts), 2, total_reads_cont_5prime, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000})  %>% 
  rownames_to_column("control_id") %>% 
  mutate(type = "5prime") %>% 
  left_join(case_controls_5prime)

names(counts_cont_5prime) <- gsub(x = names(counts_cont_5prime), pattern = "\\.bam", replacement = "")
names(counts_cont_5prime) <- gsub(x = names(counts_cont_5prime), pattern = "XX_", replacement = "control_")

counts_3prime_df <- left_join(counts_case_3prime, counts_cont_3prime)
counts_5prime_df <- left_join(counts_case_5prime, counts_cont_5prime)
counts_df <- rbind(counts_3prime_df, counts_5prime_df)
assign(paste0(day, "_counts_df"),  counts_df)
}

#Calculate and plot lfc between case and controls
d0_merge <- d0_counts_df %>% 
  pivot_longer(-c(case_id, control_id, type), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("group", "mark", "day", "replicate"), sep = "_") %>% 
  group_by(case_id, control_id, mark, day, group, type) %>% 
  dplyr::summarize(cpm = geoMean(cpm + 0.01))

d2_merge <- d2_counts_df %>% 
  pivot_longer(-c(case_id, control_id, type), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("group", "mark", "day", "replicate"), sep = "_") %>% 
  group_by(case_id, control_id, mark, day, group, type) %>% 
  dplyr::summarize(cpm = geoMean(cpm + 0.01))

d4_merge <- d4_counts_df %>% 
  pivot_longer(-c(case_id, control_id, type), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("group", "mark", "day", "replicate"), sep = "_") %>% 
  group_by(case_id, control_id, mark, day, group, type) %>% 
  dplyr::summarize(cpm = geoMean(cpm + 0.01))



counts_cnt_merge <- rbind(d0_merge, d2_merge, d4_merge) 

cnt_vec <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "H2K119Aub")

output <- data.frame()
for (m in cnt_vec) {
  for (t in c("3prime", "5prime")) {
    for (d in c("d0", "d2", "d4")) {
      a <- counts_cnt_merge %>% 
        filter(group == "case")
      
      b <- counts_cnt_merge %>% 
        filter(group == "control")
      
      p_val <- wilcox.test(a[a$mark==m&a$day==d&a$type==t,]$cpm, 
                           b[b$mark==m&b$day==d&b$type==t,]$cpm, 
                           alternative = "two.sided")$p.value
      
      row <- data.frame(mark = m, day = d, type = t, p_val = p_val)
      output <- rbind(output, row)
    }}}

pval = output

write_delim(pval,paste0(fig_dir, day, "_lfc_heatmaps_wilcoxon_pval.txt"), delim = "\t")


lfc_df <- counts_cnt_merge %>% 
  group_by(mark, day, group, type) %>% 
  summarize(cpm = geoMean(cpm)) %>% 
  pivot_wider(names_from = group, values_from = cpm) %>% 
  mutate(lfc = log2(case/control))

#Plots mean lfc per case-control pair as log2foldchange
lfc_map_3prime <- lfc_df %>%
  filter(type == "3prime") %>% 
  ggplot(aes(y = factor(mark, 
                        levels = rev(c("H3K27ac", "H3K4me3", "H3K4me1", 
                                       "H3K27me3", "H2K119Aub", "H3K9me3", "H3K36me3"))), 
             x = day)) +
  geom_tile(aes(fill = lfc), color = "black") +
  scale_fill_gradient2(low = "#0D0887", mid = "white", high = "#DE5F65", midpoint = 0, limits = c(-3, 3)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "Mean LFC (Case/Control)")

pdf(paste0(fig_dir, "Fig5D_CnT_3prime_lfc.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(lfc_map_3prime, width = unit(0.85714285714, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

lfc_map_5prime <- lfc_df %>%
  filter(type == "5prime") %>% 
  ggplot(aes(y = factor(mark, 
                        levels = rev(c("H3K27ac", "H3K4me3", "H3K4me1", 
                                       "H3K27me3", "H2K119Aub", "H3K9me3", "H3K36me3"))), 
             x = day)) +
  geom_tile(aes(fill = lfc), color = "black") +
  scale_fill_gradient2(low = "#0D0887", mid = "white", high = "#DE5F65", midpoint = 0, limits = c(-3, 3)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "Mean LFC (Case/Control)")

pdf(paste0(fig_dir, "Fig5E_CnT_5prime_lfc.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(lfc_map_5prime, width = unit(0.85714285714, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

