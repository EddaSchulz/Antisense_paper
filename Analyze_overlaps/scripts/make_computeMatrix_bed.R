#Creates control sets and groups for linegraphs of 5' overlap genes
library(tidyverse)

set.seed(22)
args <- commandArgs(trailingOnly = TRUE)
assembly_dir <-args[1] 
overlap_dir <-args[2] 
count_dir <- args[3] 
output_dir <-  args[4]
day <- args[5] 

min_length = 2000
min_ol = 1000

overlap_all <- read.delim(paste0(overlap_dir, "total_overlap_table.txt")) %>% 
  mutate(overlap_perc.y = overlap_length / length.y, overlap_perc.x = overlap_length / length.x) %>% 
  filter(overlap_length >= min_ol)

overlap.x <- overlap_all %>% 
  select(id.x, overlap_type, overlap_perc.x, length.x, ass_day = day) %>% 
  filter(between(overlap_perc.x, 0.05, 0.95) & length.x >= min_length & ass_day == day) %>% 
  transmute(id = id.x, overlap_type = overlap_type, length = length.x, overlap_perc = overlap_perc.x)

overlap.y <- overlap_all %>% 
  select(id.y, overlap_type, overlap_perc.y, length.y, ass_day = day) %>% 
  filter(between(overlap_perc.y, 0.05, 0.95) & length.y >= min_length & ass_day == day) %>% 
  transmute(id = id.y, overlap_type = overlap_type, length = length.y, overlap_perc = overlap_perc.y)

overlap_filter <- rbind(overlap.x, overlap.y) %>% 
  filter(overlap_type != "intragenic")

overlaps_5prime <- overlap_filter  %>% 
  filter(overlap_type == "5'overlap")
overlaps_3prime <- overlap_filter  %>% 
  filter(overlap_type == "3'overlap")

## Reads in all genes and filters those that don't have an overlap to act as a control
## Create BED files for plus and minus transcripts
assembly <- read.delim(paste0(assembly_dir, "XX_", day, "_trimmed_transcripts.bed"), 
                          col.names = c("chr", "start", "end", "id", "score", "strand"))

assembly_plus <- assembly %>% 
  filter(strand == "+")

assembly_minus <- assembly %>% 
  filter(strand == "-")

write_delim(assembly_plus, paste0(output_dir, "XX_", day, "_assembly_plus.bed"), delim = "\t", 
            col_names = FALSE)
write_delim(assembly_minus, paste0(output_dir, "XX_", day, "_assembly_minus.bed"), delim = "\t", 
            col_names = FALSE)

TPM_assembly <- read.delim(paste0(count_dir, "XX_", day, "_assembly_TPM_table.txt")) %>% 
  dplyr::rename(id = GeneID, strand = Strand)

free_genes <- assembly %>% 
  left_join(TPM_assembly) %>% 
  filter(!id %in% overlap_all$id.x) %>% 
  filter(!id %in% overlap_all$id.y) %>% 
  mutate(length = end - start, group = "control") %>% 
  filter(day == day & length >= min_length) %>% 
  select(-score, -chr, -start, -end)

write_delim(free_genes, paste0(output_dir, "XX_", day, "free_info.txt"), delim = "\t")

ol_list <- c("3prime", "5prime")

for (i in ol_list ){
#Prints ids of overlaps
total_case_plus <- get(paste0("overlaps_", i)) %>%
  left_join(TPM_assembly) %>% 
  mutate(group = "case")

total_case_minus <- get(paste0("overlaps_", i)) %>%
  left_join(TPM_assembly) %>% 
  mutate(group = "case")

total_case <- rbind(total_case_plus , total_case_minus) %>% 
  unique()

write_delim(total_case, paste0(output_dir, "XX_", day, "_", i, "_ids.txt"), delim = "\t")


#Calculate a set of matched controls to the 5' overlap gene
#First randomize transcripts and filter genes with very large/short overlapping regions
total_case_random <- total_case[sample(nrow(total_case)),]
loop_df <- free_genes
loop_df$TPM <- get(paste0("TPM_", day), loop_df)
case_controls <- data.frame()
for (p in 1:nrow(total_case_random)) {
  case_length = total_case_random$length[p]
  case_TPM = get(paste0("TPM_", day), total_case_random)[p]
  case_id = total_case_random$id[p]
  
  row <- loop_df %>% 
    mutate(simil_TPM = abs(log2(TPM / case_TPM)), 
           simil_length = abs(log2(length / case_length))) %>% 
    mutate(similarity = simil_TPM + simil_length) %>%
    slice_min(similarity, n = 1, with_ties = FALSE) %>%
    transmute(case_id = case_id, control_id = id, 
              control_length = length, case_length = case_length, control_TPM = TPM, case_TPM = case_TPM)
  
  loop_df <- loop_df %>% 
    filter(id != row$control_id)
  case_controls <- rbind(case_controls, row)
}


write_delim(case_controls, paste0(output_dir, "XX_", day, "_", i, "_case_controls.txt"), delim = "\t")

#Prepares BED files of complete genes
df_full <- case_controls %>% 
  pivot_longer(c(1:2), names_to = "type", values_to = "id") %>% 
  left_join(assembly)

bed_case <- df_full %>% 
  filter(type == "case_id") %>%
  select(chr, start, end, score, id, strand)

bed_cont <- df_full %>% 
  filter(type == "control_id") %>%
  select(chr, start, end, score, id, strand)

write_delim(bed_case, paste0(output_dir, "XX_", day, "_", i, "_case.bed"), delim = "\t", col_names = FALSE)
write_delim(bed_cont, paste0(output_dir, "XX_", day, "_", i, "_control.bed"), delim = "\t", col_names = FALSE)


#Prepares BED files for free regions and overlaps in cases (3')
if (i == "3prime") {
case_controls_info <- inner_join(case_controls, total_case_random, by = c("case_id" = "id")) %>% 
  left_join(assembly, by = c("case_id" = "id", "strand"))

case_free_bed_plus <- case_controls_info %>% 
  filter(strand == "+") %>% 
  mutate(overlap_start = round(end - (overlap_perc * length))) %>% 
  transmute(chr = chr, start = start, end = overlap_start, score = 1000, id = case_id, strand = strand)

case_free_bed_minus <- case_controls_info %>% 
  filter(strand == "-") %>% 
  mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
  transmute(chr = chr, start = overlap_end, end = end, score = 1000, id = case_id, strand = strand)

case_overlap_bed_plus <- case_controls_info %>% 
  filter(strand == "+") %>% 
  mutate(overlap_start = round(start + (overlap_perc * length))) %>% 
  transmute(chr = chr, start = overlap_start, end = end, score = 1000, id = case_id, strand = strand)

case_overlap_bed_minus <- case_controls_info %>% 
  filter(strand == "-") %>% 
  mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
  transmute(chr = chr, start = start, end = overlap_end, score = 1000, id = case_id, strand = strand)

#Preparing BED files for free regions and mock overlaps in controls
controls_info <- case_controls_info %>% 
  select(id = control_id, overlap_perc) %>% 
  left_join(assembly) %>% 
  mutate(length = end - start)

controls_free_bed_plus <- controls_info %>% 
  filter(strand == "+") %>% 
  mutate(overlap_start = round(end - (overlap_perc * length))) %>% 
  transmute(chr = chr, start = start, end = overlap_start, score = 1000, id = id, strand = strand)

controls_free_bed_minus <- controls_info %>% 
  filter(strand == "-") %>% 
  mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
  transmute(chr = chr, start = overlap_end, end = end, score = 1000, id = id, strand = strand)

controls_overlap_bed_plus <- controls_info %>% 
  filter(strand == "+") %>% 
  mutate(overlap_start = round(start + (overlap_perc * length))) %>% 
  transmute(chr = chr, start = overlap_start, end = end, score = 1000, id = id, strand = strand)

controls_overlap_bed_minus <- controls_info %>% 
  filter(strand == "-") %>% 
  mutate(overlap_end = round(start + (overlap_perc * length)))%>% 
  transmute(chr = chr, start = start, end = overlap_end, score = 1000, id = id, strand = strand)
} else {
  #Prepares BED files for free regions and overlaps in cases
  case_controls_info <- inner_join(case_controls, total_case_random, by = c("case_id" = "id")) %>% 
    left_join(assembly, by = c("case_id" = "id", "strand"))
  
  case_free_bed_plus <- case_controls_info %>% 
    filter(strand == "+") %>% 
    mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
    transmute(chr = chr, start = overlap_end, end = end, score = 1000, id = case_id, strand = strand)
  
  case_free_bed_minus <- case_controls_info %>% 
    filter(strand == "-") %>% 
    mutate(overlap_start = round(end - (overlap_perc * length))) %>% 
    transmute(chr = chr, start = start, end = overlap_start, score = 1000, id = case_id, strand = strand)
  
  case_overlap_bed_plus <- case_controls_info %>% 
    filter(strand == "+") %>% 
    mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
    transmute(chr = chr, start = start, end = overlap_end, score = 1000, id = case_id, strand = strand)
  
  case_overlap_bed_minus <- case_controls_info %>% 
    filter(strand == "-") %>% 
    mutate(overlap_start = round(end - (overlap_perc * length))) %>% 
    transmute(chr = chr, start = overlap_start, end = end, score = 1000, id = case_id, strand = strand)
  
  #Preparing BED files for free regions and mock overlaps in controls
  controls_info <- case_controls_info %>% 
    select(id = control_id, overlap_perc) %>% 
    left_join(assembly) %>% 
    mutate(length = end - start)
  
  controls_free_bed_plus <- controls_info %>% 
    filter(strand == "+") %>% 
    mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
    transmute(chr = chr, start = overlap_end, end = end, score = 1000, id = id, strand = strand)
  
  controls_free_bed_minus <- controls_info %>% 
    filter(strand == "-") %>% 
    mutate(overlap_start = round(end - (overlap_perc * length))) %>% 
    transmute(chr = chr, start = start, end = overlap_start, score = 1000, id = id, strand = strand)
  
  controls_overlap_bed_plus <- controls_info %>% 
    filter(strand == "+") %>% 
    mutate(overlap_end = round(start + (overlap_perc * length))) %>% 
    transmute(chr = chr, start = start, end = overlap_end, score = 1000, id = id, strand = strand)
  
  controls_overlap_bed_minus <- controls_info %>% 
    filter(strand == "-") %>% 
    mutate(overlap_start = round(end - (overlap_perc * length))) %>% 
    transmute(chr = chr, start = overlap_start, end = end, score = 1000, id = id, strand = strand)
  
}

#Printing BED files
write_delim(case_free_bed_plus, paste0(output_dir, "XX_", day, "_", i, "_case_free_plus.bed"), 
            delim = "\t", col_names = FALSE)
write_delim(case_free_bed_minus, paste0(output_dir, "XX_", day, "_", i, "_case_free_minus.bed"), 
            delim = "\t", col_names = FALSE)
write_delim(case_overlap_bed_plus, paste0(output_dir, "XX_", day, "_", i, "_case_overlap_plus.bed"), 
            delim = "\t", col_names = FALSE)
write_delim(case_overlap_bed_minus, paste0(output_dir, "XX_", day, "_", i, "_case_overlap_minus.bed"), 
            delim = "\t", col_names = FALSE)

write_delim(controls_free_bed_plus, paste0(output_dir, "XX_", day, "_", i, "_controls_free_plus.bed"), 
            delim = "\t", col_names = FALSE)
write_delim(controls_free_bed_minus, paste0(output_dir, "XX_", day, "_", i, "_controls_free_minus.bed"), 
            delim = "\t", col_names = FALSE)
write_delim(controls_overlap_bed_plus, paste0(output_dir, "XX_", day, "_", i, "_controls_overlap_plus.bed"), 
            delim = "\t", col_names = FALSE)
write_delim(controls_overlap_bed_minus, paste0(output_dir, "XX_", day, "_", i, "_controls_overlap_minus.bed"), 
            delim = "\t", col_names = FALSE)
}
