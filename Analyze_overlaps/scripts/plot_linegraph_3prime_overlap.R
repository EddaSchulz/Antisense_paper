library(tidyverse)
library(egg)
library(gridExtra)
library(viridis)
library(EnvStats)
library(stats)

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title = element_text(size = 6), 
                  axis.title.x = element_blank(),
                  axis.text = element_text(size = 6), strip.text = element_text(size = 6)))

args <- commandArgs(trailingOnly = TRUE)
data_dir = args[1]
count_dir = args[2]
fig_dir = args[3]
day = args[4]

#Plotting mock overlaps

#Read ids from bed file
case_plus_ids <- read.delim(paste0(data_dir, "XX_", day, "_3prime_case_free_plus.bed"), 
                               header = FALSE) %>% 
  select(id = V5)

case_minus_ids <- read.delim(paste0(data_dir, "XX_", day, "_3prime_case_free_minus.bed"), 
                            header = FALSE) %>% 
  select(id = V5)

controls_plus_ids <- read.delim(paste0(data_dir, "XX_", day, "_3prime_controls_free_plus.bed"), 
                            header = FALSE) %>% 
  select(id = V5)

controls_minus_ids <- read.delim(paste0(data_dir, "XX_", day, "_3prime_controls_free_minus.bed"), 
                             header = FALSE) %>% 
  select(id = V5)


#Reads all the different scaled count files for case
case_free_plus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_case_free_plus.txt"), header = FALSE, skip = 3)
colnames(case_free_plus) <- paste0("Case_", seq(1, 12))

case_free_minus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_case_free_minus.txt"), header = FALSE, skip = 3)
colnames(case_free_minus) <- paste0("Case_", seq(1, 12))

case_overlap_plus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_case_overlap_plus.txt"), header = FALSE, skip = 3)
colnames(case_overlap_plus) <- paste0("Case_", seq(13, 20))

case_overlap_minus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_case_overlap_minus.txt"), header = FALSE, skip = 3)
colnames(case_overlap_minus) <- paste0("Case_", seq(13, 20))

#Reads all the different scaled count files for controls
controls_free_plus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_controls_free_plus.txt"), header = FALSE, skip = 3)
colnames(controls_free_plus) <- paste0("controls_", seq(1, 12))

controls_free_minus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_controls_free_minus.txt"), header = FALSE, skip = 3)
colnames(controls_free_minus) <- paste0("controls_", seq(1, 12))

controls_overlap_plus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_controls_overlap_plus.txt"), header = FALSE, skip = 3)
colnames(controls_overlap_plus) <- paste0("controls_", seq(13, 20))

controls_overlap_minus <- read.delim(paste0(count_dir, "XX_", day, "_3prime_controls_overlap_minus.txt"), header = FALSE, skip = 3)
colnames(controls_overlap_minus) <- paste0("controls_", seq(13, 20))

#Binds all the count data in a single dataframe
case_counts_plus <- cbind(case_plus_ids, case_free_plus, case_overlap_plus)  %>% 
  pivot_longer(-1, names_to = "bin", values_to = "count") %>% 
  separate(bin, c("group", "bin"), sep = "_")

case_counts_minus <- cbind(case_minus_ids, case_free_minus, case_overlap_minus)  %>% 
  pivot_longer(-1, names_to = "bin", values_to = "count") %>% 
  separate(bin, c("group", "bin"), sep = "_")

controls_counts_plus <- cbind(controls_plus_ids, controls_free_plus, controls_overlap_plus)  %>% 
  pivot_longer(-1, names_to = "bin", values_to = "count") %>% 
  separate(bin, c("group", "bin"), sep = "_")

controls_counts_minus <- cbind(controls_minus_ids, controls_free_minus, controls_overlap_minus)  %>% 
  pivot_longer(-1, names_to = "bin", values_to = "count") %>% 
  separate(bin, c("group", "bin"), sep = "_")

computeMatrix_3prime <- rbind(case_counts_plus, case_counts_minus, controls_counts_plus, controls_counts_minus)

#Plots the median of the data(with significance)
order <- seq(1:20)


plot_overlap <- computeMatrix_3prime %>%
  dplyr::group_by(group, bin) %>% 
  dplyr::summarize(median = quantile(count, 0.5), low = quantile(count, 0.25), 
                   high = quantile(count, 0.75)) %>% 
  ggplot(aes(x = factor(bin, levels = order), y = median, color = group, group = group)) +
  geom_line() +
  geom_point(size = .5) +
  geom_errorbar(aes(ymin = low, ymax = high), width = .5) +
  scale_x_discrete(breaks = c(1, 13, 20), labels = c("TSS", "Overlap-Start", "TTS"), expand = c(0,3)) +
  scale_y_continuous(limits = c(0, NA), breaks = c(0, 1, 2)) +
  labs(x = "", y = "TT-seq Expression [a.u.]", color = "Group") +
  scale_color_manual(values = c("#3F2D87", "#777777"), labels = c("3'-Overlap", "No Overlaps")) +
  scale_fill_manual(values = c("#3F2D87", "#777777"), labels = c("3'-Overlap", "No Overlaps"))

pdf(paste0(fig_dir, "Fig5C_S4D_XX_", day, "_3prime_line_overlap.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(plot_overlap, width = unit(3, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

options(scipen = 100)
output <- data.frame()
for (i in 1:20) {
  a <- computeMatrix_3prime %>% 
    filter(bin == i & group == "Case")
  
  b <- computeMatrix_3prime %>% 
    filter(bin == i & group != "Case") 
  
  p_val <- wilcox.test(a$count, b$count, alternative = "two.sided")$p.value
  
  row <- data.frame(bin = i, p_val = p_val)
  output <- rbind(output, row)
}

write_delim(output, paste0(fig_dir, "XX_", day, "_line_sig_3prime_wilcoxtest.txt", delim = "\t"))
