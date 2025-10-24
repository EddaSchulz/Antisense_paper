#This script creates supplemental tables with gene name information
library(tidyverse)
library(Rgb)
library(egg)
library(gridExtra)
library(viridis)

select <- dplyr::select
day_vec <- c("d0", "d2", "d4")

args <- commandArgs(trailingOnly = TRUE)

overlap_dir  <- args[1]
files_dir <- args[2]
output_dir <- args[3]

gencode_gtf <- paste0(files_dir, "GENCODE_vM25_plus_Xert.gtf")

gene_names <- read.gtf(gencode_gtf) %>%
  select(ensembl_id = gene_id, gene = gene_name) %>%
  unique()

list_data <-lapply(day_vec, function(x)
{
  #Read data and group intragenic overlaps
  annot <-  read.delim(paste0(output_dir, "XX_", x, "_full_annotation.txt")) %>% 
    select(id, ensembl_id, gencode_type, annot) %>% 
    left_join(gene_names)
  
  annot.x <- annot %>% 
    rename(id.x = id, ensembl_id.x = ensembl_id, gencode_type.x = gencode_type, annot.x = annot, gene_name.x = gene)
  
  annot.y <- annot %>% 
    rename(id.y = id, ensembl_id.y = ensembl_id, gencode_type.y = gencode_type, annot.y = annot, gene_name.y = gene) 
  
  overlaps <- read.delim(paste0(overlap_dir, "XX_", x, "_overlap_table.txt")) %>% 
    select(id.x, id.y, overlap_type, overlap_end, overlap_start, overlap_length) %>% 
    left_join(annot.x) %>% 
    left_join(annot.y) %>% 
    mutate(day = x)
    
    write_delim(overlaps, paste0(output_dir, "XX_", x, "_annotated_overlaps.txt"), delim = "\t")
    
    return(overlaps)
  })
  
annot_overlaps <- bind_rows(list_data)

overlap_cats <- annot_overlaps %>% 
  mutate(cat = case_when(
    annot.x == "protein-coding" & annot.y == "protein-coding" ~ "both_coding",
    annot.x == "protein-coding" | annot.y == "protein-coding" ~ "one_coding",
    TRUE ~ "none_coding")) %>% 
  select(overlap_type, cat, day)

overlap_cats_percent <- overlap_cats %>%
  group_by(day, overlap_type, cat) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(day, overlap_type) %>%
  mutate(percentage = count / sum(count) * 100)

order <- c("none_coding", "one_coding", "both_coding")
type_order = c("3'overlap", "intragenic", "5'overlap")

plot_overlap_cat <- overlap_cats_percent %>% 
  ggplot() +
  facet_wrap(~factor(overlap_type, levels = type_order)) +
  geom_bar(aes(x = day, group = factor(cat, levels = order), 
               fill = factor(cat, levels = order), y = percentage), 
           color = NA, stat = "identity", position = "stack") +
  labs(y = "% of overlaps") +
  ylim(0, NA) +
  scale_fill_manual(values = c("#DEDEDE", "#676767", "#232323"))


pdf(paste0(output_dir, "FigS2G_Overlap_cat_perc.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(plot_overlap_cat, width = unit(1.5, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))
dev.off()
