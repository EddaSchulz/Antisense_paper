#Ratio area plots to see development of antisense transcription during differentiation
library(tidyverse)
library(Rsubread)
library(egg)
library(gridExtra)
library(EnvStats)
library(viridis)
library(ggsankey)


args <- commandArgs(trailingOnly = TRUE)
overlap_dir <- args[1]
bam_dir <- args[2]
fig_dir <- args[3]

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6), 
                  strip.background = element_blank(), legend.title = element_text(size = 6)))


overlaps_composite <- read.delim(paste0(overlap_dir, "XX_composite_overlap_table.txt"))


#Calculate reads within overlaps and calculate normed ratio for all timepoints
overlaps_plus_XX <- overlaps_composite %>% 
  transmute(GeneID = overlap_id, Chr = chr, Start = start, End = end, Strand = "+")

overlaps_minus_XX <- overlaps_composite %>% 
  transmute(GeneID = overlap_id, Chr = chr, Start = start, End = end, Strand = "-")

setwd(bam_dir)
temp = normalizePath(list.files(pattern="TT.*.bam$", full.names = TRUE))

overlaps_plus_XX_featureCounts <- featureCounts(temp,
                                                annot.ext = overlaps_plus_XX, isPairedEnd = TRUE, strandSpecific = 2,
                                                allowMultiOverlap = TRUE)

overlaps_minus_XX_featureCounts <- featureCounts(temp,
                                                 annot.ext = overlaps_minus_XX, isPairedEnd = TRUE, strandSpecific = 2,
                                                 allowMultiOverlap = TRUE)


total_reads <- overlaps_plus_XX_featureCounts$stat %>% 
  select(-Status) %>%
  summarise_all(funs(sum)) %>% 
  unlist()

counts_overlaps_plus <- sweep(as.matrix(overlaps_plus_XX_featureCounts$counts), 2, total_reads, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000}) %>% 
  transmute(d0_plus_r1 = TT_XX_d0_r1.bam, d0_plus_r2 = TT_XX_d0_r2.bam,
            d2_plus_r1 = TT_XX_d2_r1.bam, d2_plus_r2 = TT_XX_d2_r2.bam,
            d4_plus_r1 = TT_XX_d4_r1.bam, d4_plus_r2 = TT_XX_d4_r2.bam) %>% 
  rownames_to_column("overlap_id")

counts_overlaps_minus <- sweep(as.matrix(overlaps_minus_XX_featureCounts$counts), 2, total_reads, `/`) %>% 
  data.frame() %>% 
  mutate_all(.funs = function(x) {x * 1000000}) %>% 
  transmute(d0_minus_r1 = TT_XX_d0_r1.bam, d0_minus_r2 = TT_XX_d0_r2.bam,
            d2_minus_r1 = TT_XX_d2_r1.bam, d2_minus_r2 = TT_XX_d2_r2.bam,
            d4_minus_r1 = TT_XX_d4_r1.bam, d4_minus_r2 = TT_XX_d4_r2.bam) %>% 
  rownames_to_column("overlap_id")


counts_overlaps <- left_join(counts_overlaps_plus, counts_overlaps_minus)

pval_overlaps <- counts_overlaps %>%
  rowwise() %>%
  mutate(pval_d0 = t.test(c(d0_plus_r1, d0_plus_r2),
                       c(d0_minus_r1, d0_minus_r2), var.equal = TRUE)$p.value,
         pval_d2 = t.test(c(d2_plus_r1, d2_plus_r2),
                          c(d2_minus_r1, d2_minus_r2), var.equal = TRUE)$p.value,
         pval_d4 = t.test(c(d4_plus_r1, d4_plus_r2),
                          c(d4_minus_r1, d4_minus_r2), var.equal = TRUE)$p.value) %>%
  ungroup()  %>% 
  transmute(overlap_id, fdr_d0 = p.adjust(pval_d0, method = "BH"), 
            fdr_d2 = p.adjust(pval_d2, method = "BH"), 
            fdr_d4 = p.adjust(pval_d4, method = "BH"))



ratio_overlaps <- counts_overlaps  %>% 
  transmute(overlap_id, d0_plus = (d0_plus_r1 + d0_plus_r2) / 2, d0_minus = (d0_minus_r1 + d0_minus_r2) / 2,
            d2_plus = (d2_plus_r1 + d2_plus_r2) / 2, d2_minus = (d2_minus_r1 + d2_minus_r2) / 2,
            d4_plus = (d4_plus_r1 + d4_plus_r2) / 2, d4_minus = (d4_minus_r1 + d4_minus_r2) / 2) %>% 
  pivot_longer(c(2:7), names_to = "sample", values_to = "CPM") %>% 
  separate(sample, c("day", "strand"), sep = "_") %>% 
  pivot_wider(names_from = strand, values_from = CPM) %>% 
  mutate(plus = plus + 0.001, minus = minus + 0.001) %>% 
  mutate(overlap_ratio = (plus / (plus + minus))) %>% 
  select(-plus, -minus) %>% 
  pivot_wider(names_from = day, values_from = overlap_ratio) %>% 
  left_join(overlaps_composite) %>% 
  select(overlap_id, d0, d2, d4, overlap_type) %>% 
  left_join(pval_overlaps)



#Plotting the ratio of the composite overlaps
order <- c("3prime", "intragenic", "5prime")

ratio_day <- ratio_overlaps %>% 
  select(-fdr_d0, -fdr_d2, -fdr_d4) %>% 
  pivot_longer(-c(overlap_id, overlap_type), values_to = "ratio", names_to = "day") %>% 
  ggplot() +
  facet_wrap(~factor(overlap_type, levels = order)) +
  geom_violin(aes(x = day, y = ratio, fill = factor(overlap_type, levels = order)), color = "white") +
  geom_boxplot(aes(x = day, y = ratio, fill = factor(overlap_type, levels = order)),
               outlier.shape = NA, fill = NA, width = 0.1, color = "white", lwd = 0.1, coef = 0) + 
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.8)+
  labs(y = "Overlap Ratio [Plus / Total]")

pdf(paste0(fig_dir, "Fig4F_composite_overlap_strand_ratios.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(ratio_day, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()


#Identify switching pairs
switch_df <- ratio_overlaps %>% 
  mutate(d0_type = ifelse(d0 <= 0.35 & fdr_d0 <= 0.1, "minus_bias", 
                          ifelse(d0 >= 0.65 & fdr_d0 <= 0.1, "plus_bias", "bal")),
         d2_type = ifelse(d2 <= 0.35 & fdr_d2 <= 0.1, "minus_bias", 
                          ifelse(d2 >= 0.65 & fdr_d2 <= 0.1, "plus_bias", "bal")),
         d4_type = ifelse(d4 <= 0.35 & fdr_d4 <= 0.1, "minus_bias", 
                          ifelse(d4 >= 0.65 & fdr_d4 <= 0.1, "plus_bias", "bal")))  %>% 
  mutate(d2_switch = ifelse(d0_type == d2_type, "no_change", 
                            ifelse(d0_type == "bal" | d2_type == "bal", "change", "switch")),
         d4_switch = ifelse(d0_type == d4_type, "no_change", 
                            ifelse(d0_type == "bal" | d4_type == "bal", "change", "switch")))

switch_overlaps <- switch_df %>% 
  mutate(abs_change = abs(d4 - d0)) 
  

only_switch <- switch_overlaps %>% 
  filter(d4_switch == "switch")

write_delim(only_switch, paste0(overlap_dir, "switched_overlaps.txt"), delim = "\t")

total_overlaps_num <- switch_overlaps %>% 
  dplyr::group_by(overlap_type, d4_switch) %>% 
  dplyr::summarize(n = n()) %>% 
  pivot_wider(names_from = d4_switch, values_from = n)

write_delim(total_overlaps_num, paste0(overlap_dir, "switched_overlaps_n_d4.txt"), delim = "\t")

total_overlaps_num <- switch_overlaps %>% 
  dplyr::group_by(overlap_type, d2_switch) %>% 
  dplyr::summarize(n = n()) %>% 
  pivot_wider(names_from = d2_switch, values_from = n)

write_delim(total_overlaps_num, paste0(overlap_dir, "switched_overlaps_n_d2.txt"), delim = "\t")

#Create sankey diagrams for Day 4 
order <- c("minus_bias", "bal", "plus_bias")
ol_types <- c("intragenic", "3prime", "5prime")

for (ol in  ol_types)
{
  sankey_df <- switch_df %>%
    filter(overlap_type == ol) %>%
    make_long(d0_type, d4_type)

  sankey <- sankey_df %>% 
    ggplot(aes(x = x, 
                 next_x = next_x, 
                 node = factor(node, levels = order),  
                 next_node = factor(next_node, levels = order),
                 label = node)) +
    geom_sankey(node.color = "black", node.fill = "#DEDEDE", flow.alpha = 0.75, width = 0.4, flow.fill = "pink",
              flow.color = "black") +
    geom_sankey_text(size = 6/2.8) +
    theme_sankey(base_size = 6)

  pdf(paste0(fig_dir, "Fig_4GIK_d4_composite_overlaps_sankey_", ol, ".pdf"), useDingbats = FALSE, onefile = FALSE)

  fix <- set_panel_size(sankey, width = unit(3, "cm"), height = unit(4, "cm"))

  print(grid.arrange(fix))

  dev.off()
}

sankey_info <- switch_df %>%
  select(overlap_id, overlap_type, d0_type, d2_type, d4_type) %>% 
  pivot_longer(-c(overlap_id, overlap_type), names_to = "day", values_to = "type") %>% 
  group_by(overlap_type, day, type) %>% 
  summarize(n = n())

write_delim(sankey_info, paste0(fig_dir, "sankey_info.txt"), delim = "\t")



#Plot dot plots depending on overlap ratio
ratio_density_area_d4 <- switch_df %>% 
  ggplot() +
  facet_wrap(~overlap_type) +
  geom_point(aes(y = d0, x = d4, color = d4_switch), size = 0.2) +
  scale_color_manual(values = c("#A9A3C1", "#DEDEDE", "#3F2D87")) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.02, 0.02), limits = c(0, 1))



pdf(paste0(fig_dir, "Fig4GIK_d4_overlap_ratio_area.pdf"), useDingbats = FALSE, onefile = FALSE)

fix <- set_panel_size(ratio_density_area_d4, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()


ratio_density_area_d2 <- switch_df %>% 
  ggplot() +
  facet_wrap(~overlap_type) +
  geom_point(aes(y = d0, x = d2, color = d2_switch), size = 0.2) +
  scale_color_manual(values = c("#A9A3C1", "#DEDEDE", "#3F2D87")) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.02, 0.02), limits = c(0, 1))



pdf(paste0(fig_dir, "FigS3D-F_d2_overlap_ratio_area.pdf"), useDingbats = FALSE, onefile = FALSE)

fix <- set_panel_size(ratio_density_area_d2, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

