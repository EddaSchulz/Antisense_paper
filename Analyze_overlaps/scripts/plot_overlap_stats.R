#Plots general overlap stats
library(tidyverse)
library(egg)
library(gridExtra)
library(EnvStats)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
overlaps_dir <- args[1]
count_dir <- args[2]
fig_dir <- args[3]

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), 
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6), 
                  axis.title.x = element_blank(),
                  strip.background = element_blank(), legend.title = element_text(size = 6)))

#Read overlap files
overlaps_XX_d0 <- read.delim(paste0(overlaps_dir, "XX_d0_overlap_table.txt")) %>% 
  mutate(day = "d0")
overlaps_XX_d2 <- read.delim(paste0(overlaps_dir, "XX_d2_overlap_table.txt")) %>% 
  mutate(day = "d2")
overlaps_XX_d4 <- read.delim(paste0(overlaps_dir, "XX_d4_overlap_table.txt")) %>% 
  mutate(day = "d4")

overlaps_XX_complete <- rbind(overlaps_XX_d0, overlaps_XX_d2, overlaps_XX_d4)

#Read count files for overlap CPM
CPM_overlaps_XX_d0 <- read.delim(paste0(count_dir, "XX_d0_overlap_CPM_table.txt"))
CPM_overlaps_XX_d2 <- read.delim(paste0(count_dir, "XX_d2_overlap_CPM_table.txt"))
CPM_overlaps_XX_d4 <- read.delim(paste0(count_dir, "XX_d4_overlap_CPM_table.txt"))
  
CPM_overlaps_XX <- rbind(CPM_overlaps_XX_d0, CPM_overlaps_XX_d2, CPM_overlaps_XX_d4) 

#Join CPM with overlap data and write file to txt
overlap_all <- left_join(CPM_overlaps_XX, overlaps_XX_complete)

write.table(overlap_all, paste0(overlaps_dir, "total_overlap_table.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Create some generally descriptive plots
ratio_help <- overlap_all %>% 
  filter(overlap_length >= 500) %>%
  group_by(day, overlap_type) %>% 
  dplyr::summarize(count = n())

order <- c("3'overlap", "intragenic", "5'overlap")

ratio_day <- overlap_all %>% 
  mutate(overlap_ratio = ifelse(str_detect(overlap_id, "d0"), overlap_ratio_d0, 
                                ifelse(str_detect(overlap_id, "d2"), overlap_ratio_d2, overlap_ratio_d4))) %>%  
  filter(overlap_length >= 500) %>% 
  ggplot() +
  facet_wrap(~factor(overlap_type, levels = order)) +
  geom_violin(aes(x = day, y = overlap_ratio, fill = factor(overlap_type, levels = order)), color = "white") +
  geom_boxplot(aes(x = day, y = overlap_ratio, fill = factor(overlap_type, levels = order)),
               outlier.shape = NA, fill = NA, width = 0.1, color = "white", lwd = 0.1, coef = 0) + 
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.8)+
  labs(y = "Overlap Ratio [Plus / Total]")

pdf(paste0(fig_dir, "FigS1F_overlap_strand_ratios.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(ratio_day, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()


overlap_num_day <- overlap_all %>% 
  filter(overlap_length >= 500) %>% 
  ggplot(aes(x = day, fill = factor(overlap_type, levels = order))) +
  geom_bar(color = NA) +
  geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = .5), size = 6/2.8,
            color = "white") +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.8) +
  labs(y = "Number of Overlaps")

pdf(paste0(fig_dir, "Fig3D_overlap_num.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(overlap_num_day, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()


overlap_length_day <- overlap_all %>% 
  filter(overlap_length >= 500) %>% 
  ggplot() +
  facet_wrap(~factor(overlap_type, levels = order)) +
  geom_violin(aes(x = day, y = log10(overlap_length), fill = factor(overlap_type, levels = order)), color = "white") +
  geom_boxplot(aes(x = day, y = log10(overlap_length), fill = factor(overlap_type, levels = order)),
               outlier.shape = NA, coef = 0, fill = NA, width = 0.1, color = "white", lwd = 0.1) +
  labs(y = "Length in bp (log10)") +
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.8) +
  stat_summary(aes(y = log10(overlap_length), x = day, label = round(10^..y..)), 
                   fun=mean, geom="text", size=6/2.8) +
  scale_y_continuous(breaks = c(3,4,5))

pdf(paste0(fig_dir, "FigS1E_overlap_length.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(overlap_length_day, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

