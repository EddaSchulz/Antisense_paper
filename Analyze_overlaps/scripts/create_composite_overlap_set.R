#This script creates a composite set of 3'- and 5'overlaps between the three timepoints
library(tidyverse)
library(Rsubread)
library(bedr)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
assembly_dir <- args[1]
overlap_dir <- args[2]
bam_dir <- args[3]

cores=5

assembly_d0 <- read.delim(paste0(assembly_dir, "XX_d0_trimmed_transcripts.bed"), 
                       col.names = c("chr", "start", "end", "id", "score", "strand"))
assembly_d2 <- read.delim(paste0(assembly_dir, "XX_d2_trimmed_transcripts.bed"), 
                          col.names = c("chr", "start", "end", "id", "score", "strand"))
assembly_d4 <- read.delim(paste0(assembly_dir, "XX_d4_trimmed_transcripts.bed"), 
                          col.names = c("chr", "start", "end", "id", "score", "strand"))

assembly_total <- rbind(assembly_d0, assembly_d2, assembly_d4)

#Divide into plus and minus strand to find antisense overlaps
assembly_plus <- assembly_total %>% 
  filter(strand == "+")

assembly_minus <- assembly_total %>% 
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
         free_width.y = length.y - overlap_length)  %>% 
  filter(overlap_length >= 500)


#Quantifies reads in candidate overlaps
overlap_saf_plus <- hits_overlap  %>% 
  transmute(GeneID = overlap_id, Chr = chr.x, Start = as.numeric(overlap_start), End = as.numeric(overlap_end), 
            Strand = "+") 

overlap_saf_minus <- hits_overlap  %>% 
  transmute(GeneID = overlap_id, Chr = chr.x, Start = as.numeric(overlap_start), End = as.numeric(overlap_end), 
            Strand = "-") 


overlap_featureCounts_plus <- featureCounts(c(paste0(bam_dir, "TT_XX_d0_r1.bam"),
                                         paste0(bam_dir, "TT_XX_d0_r2.bam"),
                                         paste0(bam_dir, "TT_XX_d2_r1.bam"),
                                         paste0(bam_dir, "TT_XX_d2_r2.bam"),
                                         paste0(bam_dir, "TT_XX_d4_r1.bam"),
                                         paste0(bam_dir, "TT_XX_d4_r2.bam")),
                                        annot.ext = overlap_saf_plus, isPairedEnd = TRUE, strandSpecific = 2,
                                        nthreads = cores, allowMultiOverlap = TRUE)


overlap_featureCounts_minus <- featureCounts(c(paste0(bam_dir, "TT_XX_d0_r1.bam"),
                                              paste0(bam_dir, "TT_XX_d0_r2.bam"),
                                              paste0(bam_dir, "TT_XX_d2_r1.bam"),
                                              paste0(bam_dir, "TT_XX_d2_r2.bam"),
                                              paste0(bam_dir, "TT_XX_d4_r1.bam"),
                                              paste0(bam_dir, "TT_XX_d4_r2.bam")),
                                            annot.ext = overlap_saf_minus, isPairedEnd = TRUE, strandSpecific = 2,
                                            nthreads = cores, allowMultiOverlap = TRUE)

total_reads <- overlap_featureCounts_plus$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_overlaps_plus <- sweep(as.matrix(overlap_featureCounts_plus$counts), 2, total_reads, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000}) %>% 
  transmute(d0_plus = (TT_XX_d0_r1.bam + TT_XX_d0_r2.bam / 2), d2_plus = (TT_XX_d2_r1.bam + TT_XX_d2_r2.bam / 2),
            d4_plus = (TT_XX_d4_r1.bam + TT_XX_d4_r2.bam / 2)) %>% 
  rownames_to_column("overlap_id")

counts_overlaps_minus <- sweep(as.matrix(overlap_featureCounts_minus$counts), 2, total_reads, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000}) %>% 
  transmute(d0_minus = (TT_XX_d0_r1.bam + TT_XX_d0_r2.bam / 2), d2_minus = (TT_XX_d2_r1.bam + TT_XX_d2_r2.bam / 2),
            d4_minus = (TT_XX_d4_r1.bam + TT_XX_d4_r2.bam / 2)) %>% 
  rownames_to_column("overlap_id")


counts_overlaps <- left_join(counts_overlaps_plus, counts_overlaps_minus)

#Filters based on the overlap having expression in all timepoints (from either strand)
counts_overlaps_filter <- counts_overlaps %>% 
  filter(d0_plus + d0_minus >= 3 & d2_plus + d2_minus >= 3 & d4_plus + d4_minus >= 3)

hits_filtered <- left_join(counts_overlaps_filter, hits_overlap) 

#Splits dfs according to overlap type
hits_3prime <- hits_filtered %>% 
  filter(overlap_type == "3'overlap") %>% 
  select(chr = chr.x, start = overlap_start, end = overlap_end)


hits_5prime <- hits_filtered %>% 
  filter(overlap_type == "5'overlap") %>% 
  select(chr = chr.x, start = overlap_start, end = overlap_end)


hits_intragenic <- hits_filtered %>% 
  filter(overlap_type == "intragenic") %>% 
  select(chr = chr.x, start = overlap_start, end = overlap_end)


#Merge overlaps together

bedr_5prime <- hits_5prime  %>% 
  transmute(region = paste0(chr, ":", as.numeric(start), "-", as.numeric(end))) %>% 
  unlist(use.names = FALSE)

bedr_5prime_sorted <- bedr.sort.region(bedr_5prime, method = "natural")
bedr_5prime_merged <- bedr.merge.region(bedr_5prime_sorted)
bedr_df_5prime <- data.frame(region = bedr_5prime_merged) %>% 
  mutate(overlap_type = "5prime")


bedr_intragenic <- hits_intragenic  %>% 
  transmute(region = paste0(chr, ":", as.numeric(start), "-", as.numeric(end))) %>% 
  unlist(use.names = FALSE)

bedr_intragenic_sorted <- bedr.sort.region(bedr_intragenic, method = "natural")
bedr_intragenic_merged <- bedr.merge.region(bedr_intragenic_sorted)
bedr_intragenic_filter_5prime <- bedr.join.region(bedr_intragenic_merged, bedr_5prime_merged) %>% 
  filter(V4 != ".")
bedr_df_intragenic <- data.frame(region = bedr_intragenic_merged) %>% 
  mutate(overlap_type = "intragenic") %>% 
  filter(!region %in% bedr_intragenic_filter_5prime$index)


bedr_3prime <- hits_3prime  %>% 
  transmute(region = paste0(chr, ":", as.numeric(start), "-", as.numeric(end))) %>% 
  unlist(use.names = FALSE)

bedr_3prime_sorted <- bedr.sort.region(bedr_3prime, method = "natural")
bedr_3prime_merged <- bedr.merge.region(bedr_3prime_sorted)
bedr_3prime_filter_5prime <- bedr.join.region(bedr_3prime_merged, bedr_5prime_merged) %>% 
  filter(V4 != ".")
bedr_3prime_filter_intragenic <- bedr.join.region(bedr_3prime_merged, bedr_intragenic_merged) %>% 
  filter(V4 != ".")  
bedr_df_3prime <- data.frame(region = bedr_3prime_merged) %>% 
  mutate(overlap_type = "3prime") %>% 
  filter(!region %in% bedr_3prime_filter_5prime$index) %>% 
  filter(!region %in% bedr_3prime_filter_intragenic$index)




#Give new IDs and print to file
bedr_df <- rbind(bedr_df_3prime, bedr_df_intragenic, bedr_df_5prime) %>% 
  mutate(overlap_id = paste0("COMP_OVERLAP_", region)) %>% 
  separate(region, c("chr", "region"), sep = ":") %>% 
  separate(region, c("start", "end"), sep = "-") %>% 
  mutate(overlap_length = as.numeric(end) - as.numeric(start)) %>% 
  select(overlap_id, everything()) %>% 
  filter(overlap_length >= 500)


write_delim(bedr_df, paste0(overlap_dir, "XX_composite_overlap_table.txt"), delim = "\t")

bedr_bed <- bedr_df %>% 
  transmute(chr, start, end, overlap_id, score = 1000, strand = ".") 
  
            
write_delim(bedr_bed, paste0(overlap_dir, "XX_composite_overlap_table.bed"), delim = "\t", col_names = FALSE)



