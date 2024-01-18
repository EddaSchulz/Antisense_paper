library(tidyverse)
chr_vec <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

args <- commandArgs(trailingOnly = TRUE)
day <- args[1] 
input_dir <- args[2]
output_dir <- args[3]

#TSS

binned_counts_tss_plus <- read.delim(paste0(input_dir, "XX_", day ,"_tss_plus_binned_counts.txt")) %>% 
  arrange(factor(chr, levels = chr_vec), start)

plus_transcripts_list <- read.delim(paste0(input_dir, "XX_", day ,"_plus_transcripts_list.bed"), 
                                    header = FALSE, col.names = c("chr", "start", "end", "gene", "gene_type",
                                                                  "trans_nr")) %>% 
  arrange(factor(chr, levels = chr_vec), start) %>% 
  mutate(range = (end - start) / 40, -1)

trans_vec_plus <- unlist(mapply(function(x, y) {rep(x, each = y)}, 
                         plus_transcripts_list$trans_nr, plus_transcripts_list$range))

tss_plus_starts <- binned_counts_tss_plus %>% 
  cbind(trans_vec_plus)

tss_plus_help <- tss_plus_starts %>% 
  group_by(trans_vec_plus) %>% 
  dplyr::summarize(cut = mean(TT_plus) / 4, max = max(TT_plus))
  
tss_plus_starts2 <- left_join(tss_plus_starts, tss_plus_help) %>% 
  filter(TT_plus >= cut) %>%
  group_by(trans_vec_plus) %>% 
  slice_min(start)  %>% 
  arrange(factor(chr, levels = chr_vec), start) %>% 
  select(chr, start, trans_vec_plus)

#TTS

binned_counts_tts_plus <- read.delim(paste0(input_dir, "XX_", day ,"_tts_plus_binned_counts.txt")) %>% 
  arrange(factor(chr, levels = chr_vec), start)

tts_plus_ends <- binned_counts_tts_plus %>% 
  cbind(trans_vec_plus)

tts_plus_help <- tts_plus_ends %>% 
  group_by(trans_vec_plus) %>% 
  dplyr::summarize(cut = mean(TT_plus) / 2, max = max(TT_plus))

tts_plus_ends2 <- left_join(tts_plus_ends, tts_plus_help) %>% 
  filter(TT_plus >= cut) %>%
  group_by(trans_vec_plus) %>% 
  slice_max(end)  %>% 
  arrange(factor(chr, levels = chr_vec), start) %>% 
  select(chr, end, trans_vec_plus)

#Combine
trimmed_plus <- left_join(tss_plus_starts2, tts_plus_ends2) %>% 
  select(chr, start, end, id = trans_vec_plus)

##MINUS

#TTS

binned_counts_tts_minus <- read.delim(paste0(input_dir, "XX_", day ,"_tts_minus_binned_counts.txt")) %>% 
  arrange(factor(chr, levels = chr_vec), start)

minus_transcripts_list <- read.delim(paste0(input_dir, "XX_", day ,"_minus_transcripts_list.bed"), 
                                    header = FALSE, col.names = c("chr", "start", "end", "gene", "gene_type",
                                                                  "trans_nr"))  %>% 
  arrange(factor(chr, levels = chr_vec), start) %>% 
  mutate(range = (end - start) / 40)

trans_vec_minus <- unlist(mapply(function(x, y) {rep(x, each = y)}, 
                                minus_transcripts_list$trans_nr, minus_transcripts_list$range))

tts_minus_starts <- binned_counts_tts_minus %>% 
  cbind(trans_vec_minus)

tts_minus_help <- tts_minus_starts %>% 
  group_by(trans_vec_minus) %>% 
  dplyr::summarize(cut = mean(TT_minus) / 2, max = max(TT_minus))

tts_minus_starts2 <- left_join(tts_minus_starts, tts_minus_help) %>% 
  filter(TT_minus >= cut) %>%
  group_by(trans_vec_minus) %>% 
  slice_min(start)  %>% 
  arrange(factor(chr, levels = chr_vec), start) %>% 
  select(chr, start, trans_vec_minus)

#TSS

binned_counts_tss_minus <- read.delim(paste0(input_dir, "XX_", day ,"_tss_minus_binned_counts.txt")) %>% 
  arrange(factor(chr, levels = chr_vec), start)

tss_minus_ends <- binned_counts_tss_minus %>% 
  cbind(trans_vec_minus)

tss_minus_help <- tss_minus_ends %>% 
  group_by(trans_vec_minus) %>% 
  dplyr::summarize(cut = mean(TT_minus) / 4, max = max(TT_minus))

tss_minus_ends2 <- left_join(tss_minus_ends, tss_minus_help) %>% 
  filter(TT_minus >= cut) %>%
  group_by(trans_vec_minus) %>% 
  slice_max(end)  %>% 
  arrange(factor(chr, levels = chr_vec), start) %>% 
  select(chr, end, trans_vec_minus)

#Combine
trimmed_minus <- left_join(tts_minus_starts2, tss_minus_ends2) %>% 
  select(chr, start, end, id = trans_vec_minus)

#Combine plus and minus
trimmed_total <- rbind(trimmed_plus, trimmed_minus) %>% 
  mutate(score = 0, strand = ifelse(str_detect(id, "plus"), "+", "-")) %>% 
  arrange(factor(chr, levels = chr_vec), start)

write.table(trimmed_total, file=paste0(output_dir, "XX_", day,"_trimmed_transcripts.bed"), quote=F, 
            sep="\t", row.names=F, col.names=F)

