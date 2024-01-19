#This script counts TT-seq reads in the complete annotation and within the overlaps
library(tidyverse)
library(Rsubread)
library(egg)
library(gridExtra)
library(stats)


args <- commandArgs(trailingOnly = TRUE)
day <- args[1]
bam_dir <-  args[2]
overlaps_dir <- args[3]
assembly_dir <- args[4]
output_dir <- args[5]

assembly <- read.delim(paste0(assembly_dir, "XX_", day, "_trimmed_transcripts.bed"), 
                       col.names = c("chr", "start", "end", "id", "score", "strand"))

overlaps <- read.delim(paste0(overlaps_dir, "XX_", day, "_overlap_table.txt"))

assembly_saf <- assembly %>% 
  select(GeneID = id, Chr = chr, Start = start, End = end, Strand = strand)


assembly_featureCounts <- featureCounts(c(paste0(bam_dir, "TT_XX_d0_r1.bam"),
                                        paste0(bam_dir, "TT_XX_d0_r2.bam"),
                                        paste0(bam_dir, "TT_XX_d2_r1.bam"),
                                        paste0(bam_dir, "TT_XX_d2_r2.bam"),
                                        paste0(bam_dir, "TT_XX_d4_r1.bam"),
                                        paste0(bam_dir, "TT_XX_d4_r2.bam")),
                                        annot.ext = assembly_saf, isPairedEnd = TRUE, strandSpecific = 2)

counts_assembly <- data.frame(assembly_featureCounts$counts, assembly_featureCounts$annotation)

names(counts_assembly) <- gsub(x = names(counts_assembly), pattern = "\\.", replacement = "_")
names(counts_assembly) <- gsub(x = names(counts_assembly), pattern = "TT_", replacement = "")
names(counts_assembly) <- gsub(x = names(counts_assembly), pattern = "_bam", replacement = "")

TPM_assembly <- counts_assembly %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length, -Chr, -Start, -End, -Strand), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length)

TPM_assembly$TPM_d0 = (unlist(TPM_assembly[,1]) + unlist(TPM_assembly[,2])) / 2
TPM_assembly$TPM_d2 = (unlist(TPM_assembly[,3]) + unlist(TPM_assembly[,4])) / 2
TPM_assembly$TPM_d4 = (unlist(TPM_assembly[,5]) + unlist(TPM_assembly[,6])) / 2

#Calculate reads within overlaps
overlaps_plus <- overlaps %>% 
  transmute(GeneID = overlap_id, Chr = chr.x, Start = overlap_start, End = overlap_end, Strand = "+")

overlaps_minus <- overlaps %>% 
  transmute(GeneID = overlap_id, Chr = chr.x, Start = overlap_start, End = overlap_end, Strand = "-")

overlaps_plus_featureCounts <- featureCounts(c(paste0(bam_dir, "TT_XX_d0_r1.bam"),
                                               paste0(bam_dir, "TT_XX_d0_r2.bam"),
                                               paste0(bam_dir, "TT_XX_d2_r1.bam"),
                                               paste0(bam_dir, "TT_XX_d2_r2.bam"),
                                               paste0(bam_dir, "TT_XX_d4_r1.bam"),
                                               paste0(bam_dir, "TT_XX_d4_r2.bam")),
                                        annot.ext = overlaps_plus, isPairedEnd = TRUE, strandSpecific = 2,
                                        allowMultiOverlap = TRUE)

overlaps_minus_featureCounts <- featureCounts(c(paste0(bam_dir, "TT_XX_d0_r1.bam"),
                                                paste0(bam_dir, "TT_XX_d0_r2.bam"),
                                                paste0(bam_dir, "TT_XX_d2_r1.bam"),
                                                paste0(bam_dir, "TT_XX_d2_r2.bam"),
                                                paste0(bam_dir, "TT_XX_d4_r1.bam"),
                                                paste0(bam_dir, "TT_XX_d4_r2.bam")),
                                             annot.ext = overlaps_minus, isPairedEnd = TRUE, strandSpecific = 2,
                                             allowMultiOverlap = TRUE)

total_reads <- overlaps_plus_featureCounts$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_overlaps_plus <- data.frame(sweep(as.matrix(overlaps_plus_featureCounts$counts), 2, total_reads, `/`)) %>%
  mutate_all(.funs = function(x) {x * 1000000}) %>% 
  transmute(d0_plus = (TT_XX_d0_r1.bam + TT_XX_d0_r2.bam / 2), d2_plus = (TT_XX_d2_r1.bam + TT_XX_d2_r2.bam / 2),
            d4_plus = (TT_XX_d4_r1.bam + TT_XX_d4_r2.bam / 2)) %>% 
  rownames_to_column("overlap_id")

counts_overlaps_minus <- data.frame(sweep(as.matrix(overlaps_minus_featureCounts$counts), 2, total_reads, `/`)) %>% 
  mutate_all(.funs = function(x) {x * 1000000}) %>% 
  transmute(d0_minus = (TT_XX_d0_r1.bam + TT_XX_d0_r2.bam / 2), d2_minus = (TT_XX_d2_r1.bam + TT_XX_d2_r2.bam / 2),
            d4_minus = (TT_XX_d4_r1.bam + TT_XX_d4_r2.bam / 2)) %>% 
  rownames_to_column("overlap_id")


counts_overlaps <- left_join(counts_overlaps_plus, counts_overlaps_minus)


#Add counts to large table
Table_assembly <- TPM_assembly %>% 
  select(GeneID, Strand, TPM_d0, TPM_d2, TPM_d4)

Table_overlaps <- counts_overlaps %>% 
  mutate(overlap_ratio_d0 = (d0_plus / (d0_minus + d0_plus)),
         overlap_ratio_d2 = (d2_plus / (d2_minus + d2_plus)),
         overlap_ratio_d4 = (d4_plus / (d4_minus + d4_plus)))
  

write_delim(Table_assembly, paste0(output_dir, "XX_", day, "_assembly_TPM_table.txt"), delim = "\t")
write_delim(Table_overlaps, paste0(output_dir, "XX_", day, "_overlap_CPM_table.txt"), delim = "\t")



