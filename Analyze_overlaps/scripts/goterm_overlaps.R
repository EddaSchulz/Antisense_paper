#This script does GO Term enrichment for the different overlap groups
library(tidyverse)
library(egg)
library(gridExtra)
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)

select <- dplyr::select
day_vec <- c("d0", "d2", "d4")

args <- commandArgs(trailingOnly = TRUE)
overlap_dir  <- args[1]
output_dir <- args[2]


list_out <- lapply(day_vec, function(x)
{
  #Read data and group intragenic overlaps
  annot <-  read.delim(paste0(output_dir, "XX_", x, "_full_annotation.txt"))
  overlaps <- read.delim(paste0(overlap_dir, "XX_", x, "_overlap_table.txt")) %>% 
    select(id.x, id.y, overlap_type) %>% 
    pivot_longer(c(id.x, id.y), values_to = "id", names_to = "strand") %>% 
    select(-strand, -free_width.x, -free_width.y) %>% 
    unique() %>% 
    left_join(annot) %>% 
    mutate(ensembl_id = str_extract(ensembl_id, "^[^.]+"))
  
  #Get overlap GENCODE sets and group as list
  prom_overlap <- overlaps %>% 
    filter(overlap_type == "5'overlap") %>% 
    filter(GENCODE == "TRUE")
  
  end_overlap <-  overlaps %>% 
    filter(overlap_type == "3'overlap") %>% 
    filter(GENCODE == "TRUE")
  
  free <- annot %>% 
    filter(!id %in% overlaps$id) %>% 
    filter(GENCODE == "TRUE")%>% 
    mutate(ensembl_id = str_extract(ensembl_id, "^[^.]+"))
  
  #Specify background set (all detected transcripts)
  entrez_bg <- mapIds(org.Mm.eg.db, keys = str_extract(annot[annot$GENCODE=="TRUE",]$ensembl_id, "^[^.]+"),
                       column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  
  #Perform GO Term enrichment
  go_function <- function(x) {
    entrez_ids <- mapIds(org.Mm.eg.db, keys = x,
                         column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
    enrich_out <- enrichGO(gene = entrez_ids, 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "fdr",
                             universe = entrez_bg,
                            pvalueCutoff = 1)
    
    return(enrich_out)
  }
  
  
  prom_enrich <- go_function(prom_overlap$ensembl_id)
  end_enrich <- go_function(end_overlap$ensembl_id)
  free_enrich <- go_function(free$ensembl_id)
  
  #Get out results and store in dataframe
  prom_df <- prom_enrich@result %>% 
    mutate(overlap_type = "prom_overlap")
  
  end_df <- end_enrich@result %>% 
    mutate(overlap_type = "3'overlap")
  
  free_df <- free_enrich@result %>% 
    mutate(overlap_type = "free")
  
  combi_go <- rbind(prom_df, end_df, free_df) %>% 
    mutate(day = x)
  
  return(combi_go)

})

go_df <- bind_rows(list_out)

sig_terms <- go_df %>% 
  select(ID, Description, p.adjust) %>%  
  filter(p.adjust <= 0.05) %>% 
  select(-p.adjust) %>% 
  unique()

plot_df <- go_df %>% 
  filter(ID %in% sig_terms$ID)

# Convert data to a matrix for clustering
clust_matrix <- plot_df %>%
  mutate(type_day = paste0(overlap_type, "_", day)) %>% 
  select(type_day, Description, p.adjust) %>%
  pivot_wider(names_from = type_day, values_from = p.adjust) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# Perform hierarchical clustering
clust_order <- hclust(dist(clust_matrix))$order

# Reorder Description based on clustering
plot_df <- plot_df %>%
  mutate(Description = factor(Description, levels = rownames(clust_matrix)[clust_order]))

ol_order <- c("3'overlap", "prom_overlap", "free")

dotplot <- plot_df %>% 
  ggplot(aes(x = day, 
             y = Description, size = -log10(p.adjust), fill = log2(FoldEnrichment))) +
  facet_wrap(~factor(overlap_type, levels = ol_order)) +
  geom_point(shape = 21) +
  scale_size_continuous(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05), -log10(0.2)), range = c(0.5, 5)) +
  scale_fill_gradient2(low = "#0D0887", mid = "white", high = "#DE5F65", midpoint = 0, limits = c(-2, 2)) 

pdf(paste0(output_dir, "FigS3A_GO_terms_total.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(dotplot, width = unit(1.5, "cm"), height = unit(10, "cm"))

print(grid.arrange(fix))
dev.off()

