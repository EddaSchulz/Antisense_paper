#This script detects all overlaps in a TT-seq annotation
library(tidyverse)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
day <- args[1] 
input_dir <- args[2]
output_dir <- args[3]


assembly <- read.delim(paste0(input_dir, "XX_", day, "_trimmed_transcripts.bed"), 
                       col.names = c("chr", "start", "end", "id", "score", "strand"))

#Divide into plus and minus strand to find antisense overlaps
assembly_plus <- assembly %>% 
  filter(strand == "+")

assembly_minus <- assembly %>% 
  filter(strand == "-")

generanges_plus <- makeGRangesFromDataFrame(assembly_plus, keep.extra.columns = TRUE)

generanges_minus <- makeGRangesFromDataFrame(assembly_minus, keep.extra.columns = TRUE)

hits <- findOverlaps(generanges_plus, generanges_minus, ignore.strand = TRUE) 

#Puts overlap in dataframe and collects gene info for both strands
plus_hits <- as.data.frame(generanges_plus[queryHits(hits)]) %>% 
  select(chr.x = seqnames, start.x = start, end.x = end, strand.x = strand, id.x = id, length.x = width)

minus_hits <- as.data.frame(generanges_minus[subjectHits(hits)]) %>% 
  dplyr::select(chr.y = seqnames, start.y = start, end.y = end, strand.y = strand, id.y = id, length.y = width)

hits_df <- cbind(plus_hits, minus_hits)

#Finds start and end of overlapping regions
hits_df$overlap_end <- apply(hits_df, 1, function(x) sort(c(x[2], x[3], x[8], x[9]), partial = 3)[3])
hits_df$overlap_start <- apply(hits_df, 1, function(x) sort(c(x[2], x[3], x[8], x[9]), partial = 2)[2])

#Groups overlaps according to gene architecture
hits_overlap <- hits_df %>% 
  mutate(overlap_length = abs(as.numeric(overlap_start) - as.numeric(overlap_end)) + 1) %>% 
  mutate(overlap_type = ifelse(start.x < start.y & end.x < end.y, "3'overlap", 
                               ifelse(start.x > start.y & end.x > end.y, "5'overlap", "intragenic"))) %>% 
  mutate(overlap_id = paste0(id.x, ":", id.y), free_width.x = length.x - overlap_length, 
         free_width.y = length.y - overlap_length) 

write_delim(hits_overlap, paste0(output_dir, "XX_", day, "_overlap_table.txt"), delim = "\t")
