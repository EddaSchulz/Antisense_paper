#This script plots statistics for the control regions for 3' and 5'
library(tidyverse)
library(egg)
library(gridExtra)
library(stats)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
assembly_dir <- args[2]
count_dir <- args[3]
fig_dir <- args[4]
day <- args[5]

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), 
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title.x = element_blank(),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6), 
                  strip.background = element_blank(), legend.title = element_text(size = 6)))


assembly_TPM <- read.delim(paste0(count_dir, "XX_", day, "_assembly_TPM_table.txt")) %>% 
  transmute(id = GeneID, TPM_d0, TPM_d2, TPM_d4)

case_5prime <- read.delim(paste0(data_dir, "XX_", day, "_5prime_case.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand")) %>% 
  transmute(group = "5prime", type = "5prime_case", id = id, length = end - start)

control_5prime <- read.delim(paste0(data_dir, "XX_", day, "_5prime_control.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand")) %>% 
  transmute(group = "5prime", type = "5prime_control", id = id, length = end - start)

case_3prime <- read.delim(paste0(data_dir, "XX_", day, "_3prime_case.bed"), header = FALSE, 
                          col.names = c("chr", "start", "end", "score", "id", "strand")) %>% 
  transmute(group = "3prime", type = "3prime_case", id = id, length = end - start)

control_3prime <- read.delim(paste0(data_dir, "XX_", day, "_3prime_control.bed"), header = FALSE, 
                             col.names = c("chr", "start", "end", "score", "id", "strand")) %>% 
  transmute(group = "3prime", type = "3prime_control", id = id, length = end - start)

stats_df <- rbind(case_5prime, control_5prime, case_3prime, control_3prime) %>% 
  left_join(assembly_TPM) %>% 
  pivot_longer(-c(1:3), names_to = "stat", values_to = "value")

stat_violin_tpm <- stats_df %>% 
  filter(stat != "length") %>% 
  mutate(value = log2(value + 0.01)) %>% 
  ggplot(aes(x = stat, y = value, fill = type)) +
  facet_wrap(~group) +
  geom_violin(color = "white", position = position_dodge(width = 0.9)) +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1, color = "white", lwd = 0.1, 
               position = position_dodge(width = 0.9)) +
  labs(y = "TPM + 0.01 [log2]") +
  scale_fill_manual(values = c("#3F2D87", "#777777", "#F7A539", "#777777")) +
  scale_y_continuous(limits = c(-2, 10), breaks = c(0, 5, 10))

pdf(paste0(fig_dir, "Fig4B_S2C_", day, "_mock_control_expression.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(stat_violin_tpm, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()


stat_violin_len <- stats_df %>% 
  filter(stat == "length") %>% 
  mutate(value = log10(value)) %>% 
  ggplot(aes(x = stat, y = value, fill = type)) +
  facet_wrap(~group) +
  geom_violin(color = "white", position = position_dodge(width = 0.9)) +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1, color = "white", lwd = 0.1, 
               position = position_dodge(width = 0.9)) +
  labs(y = "Length in bp (log10)") +
  scale_fill_manual(values = c("#3F2D87", "#777777", "#F7A539", "#777777")) +
  scale_y_continuous(limits = c(2.5, 5.5), breaks = c(3, 4, 5))

pdf(paste0(fig_dir, "FigS2A-B_", day, "_mock_control_length.pdf"), useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(stat_violin_len, width = unit(1, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()




stats_vec <- c("TPM_d0", "TPM_d2", "TPM_d4", "length")

output <- data.frame()
for (stat in stats_vec) {
  a <- stats_df %>% 
    filter(type == "5prime_case")
  
  b <- stats_df %>% 
    filter(type == "5prime_control")
  
  p_val <- wilcox.test(a[a$stat==stat,]$value, b[b$stat==stat,]$value, alternative = "two.sided")$p.value
  
  row <- data.frame(stat = stat, p_val = p_val, comp = "5prime")
  output <- rbind(output, row)
}

pval_5prime = output

output <- data.frame()
for (stat in stats_vec) {
  a <- stats_df %>% 
    filter(type == "3prime_case")
  
  b <- stats_df %>% 
    filter(type == "3prime_control")
  
  p_val <- wilcox.test(a[a$stat==stat,]$value, b[b$stat==stat,]$value, alternative = "two.sided")$p.value
  
  row <- data.frame(stat = stat, p_val = p_val, comp = "3prime")
  output <- rbind(output, row)
}

pval_3prime = output

pval <- rbind(pval_3prime, pval_5prime)

write_delim(pval,paste0(data_dir, day, "_mock_control_comp_pval.txt"), delim = "\t")







