#Script is used to find the source of each transcript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
day <- args[1]
input_dir <- args[2]
output_dir <- args[3]


#Looks where the TSS on the plus strand were sourced from (GENCODE, ChromHMM, FANTOM5 CAGE)
trans_list_plus <- read.delim(paste0(input_dir, "XX_", day,"_plus_transcripts_list.bed"), 
                         col.names = c("chr", "start", "end", "GENCODE", "type", "id"))

chromHMM_tss_plus <- read.delim(paste0(input_dir, day,"_XX_pot_tss_plus_chromHMM_overlap.bed"), 
                           col.names = c("chr", "start", "end", "chromHMM")) %>% 
  select(-end)

CAGE_tss_plus <- read.delim(paste0(input_dir, day,"_XX_pot_tss_plus_CAGE_overlap.bed"), 
                            col.names = c("chr", "start", "end", "CAGE")) %>% 
  select(-end)

trans_sources_plus <- left_join(trans_list_plus, chromHMM_tss_plus) %>% 
  left_join(CAGE_tss_plus)

#Looks where the TSS on the minus strand were sourced from (GENCODE, ChromHMM, FANTOM5 CAGE)
trans_list_minus <- read.delim(paste0(input_dir, "XX_", day,"_minus_transcripts_list.bed"), 
                              col.names = c("chr", "start", "end", "GENCODE", "type", "id"))

chromHMM_tss_minus <- read.delim(paste0(input_dir, day,"_XX_pot_tss_minus_chromHMM_overlap.bed"), 
                                col.names = c("chr", "start", "end", "chromHMM")) %>% 
  select(-start)

CAGE_tss_minus <- read.delim(paste0(input_dir, day,"_XX_pot_tss_minus_CAGE_overlap.bed"), 
                            col.names = c("chr", "start", "end", "CAGE")) %>% 
  select(-start)

trans_sources_minus <- left_join(trans_list_minus, chromHMM_tss_minus) %>% 
  left_join(CAGE_tss_minus)


#Combines both dataframes and removes unnecessary data
trans_sources <- rbind(trans_sources_plus, trans_sources_minus) %>% 
  transmute(id, GENCODE = ifelse(GENCODE == ".", FALSE, TRUE), 
            chromHMM = ifelse(is.na(chromHMM), FALSE, as.character(chromHMM)),
            CAGE = ifelse(is.na(CAGE), FALSE, TRUE)) %>% 
  unique()

write.table(trans_sources, file=paste0(output_dir, day, "_XX_transcripts_source.txt"), quote=F, 
            sep="\t", row.names=F, col.names=T)