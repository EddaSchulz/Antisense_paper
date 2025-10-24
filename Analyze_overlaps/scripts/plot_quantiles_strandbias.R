#Describe strand bias vs total expression
library(tidyverse)
library(egg)
library(gridExtra)
library(viridis)

day_vec <- c("d0", "d2", "d4")

args <- commandArgs(trailingOnly = TRUE)
overlap_dir <- args[1]
output_dir <- args[2]

overlaps <- paste0(overlap_dir, "total_overlap_table.txt")

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), 
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6), 
                  axis.title.x = element_blank(),
                  strip.background = element_blank(), legend.title = element_text(size = 6)))

list_out <- lapply(day_vec, function(x) {
  overlap_df <- overlaps %>% 
    filter(overlap_length >= 500) %>%
    filter(day == x) %>% 
    select(overlap_id, plus_CPM = paste0(x, "_plus"), minus_CPM = paste0(x, "_minus"), 
           overlap_ratio = paste0("overlap_ratio_", x), overlap_type, day, overlap_length) %>% 
    mutate(total_CPM = plus_CPM + minus_CPM) %>% 
    mutate(total_RPKM = total_CPM / (overlap_length / 1000)) %>% 
    group_by(overlap_type) %>% 
    mutate(quantile = ntile(total_CPM, 4))
  
  return(overlap_df)
})
  

plot_df <- bind_rows(list_out)


order <- c("3'overlap", "intragenic", "5'overlap")

ratio_exp <- plot_df %>% 
  ggplot() +
  facet_wrap(day~factor(overlap_type, levels = order)) +
  geom_violin(aes(x = quantile, group = quantile,
                  y = overlap_ratio, fill = factor(overlap_type, levels = order)), color = "white") +
  geom_boxplot(aes(x = quantile, group = quantile,
                   y = overlap_ratio, fill = factor(overlap_type, levels = order)),
               outlier.shape = NA, fill = NA, width = 0.1, color = "white", lwd = 0.1, coef = 0) + 
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.8)+
  labs(y = "Overlap Ratio [Plus / Total]")

pdf(paste0(output_dir, "FigS3C_total_CPM_strandbias.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(ratio_exp, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

