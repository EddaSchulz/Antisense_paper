library(tidyverse)
library(egg)
library(gridExtra)
library(UpSetR)
library(viridis)

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), 
                  legend.text = element_text(size = 6),
                  axis.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.title.x = element_blank()))

args <- commandArgs(trailingOnly = TRUE)
files_dir <- args[1]
output_dir <- args[2]

assembly_d0 <-  read.delim(paste0(files_dir, "XX_d0_trimmed_transcripts.bed"), 
                           col.names = c("chr", "start", "end", "id", "score", "strand"))

assembly_d2 <-  read.delim(paste0(files_dir, "XX_d2_trimmed_transcripts.bed"), 
                           col.names = c("chr", "start", "end", "id", "score", "strand"))

assembly_d4 <-  read.delim(paste0(files_dir, "XX_d4_trimmed_transcripts.bed"), 
                           col.names = c("chr", "start", "end", "id", "score", "strand"))

source_d0 <- read.delim(paste0(files_dir, "d0_XX_transcripts_source.txt")) %>% 
  filter(id %in% assembly_d0$id) %>% 
  mutate(day = "d0")

source_d2 <- read.delim(paste0(files_dir, "d2_XX_transcripts_source.txt")) %>% 
  filter(id %in% assembly_d2$id) %>% 
  mutate(day = "d2")

source_d4 <- read.delim(paste0(files_dir, "d4_XX_transcripts_source.txt")) %>% 
  filter(id %in% assembly_d4$id) %>% 
  mutate(day = "d4") 

total_sources <- rbind(source_d0, source_d2, source_d4) %>% 
  mutate(chromHMM_prom = ifelse(chromHMM == "promoter", TRUE, FALSE),
         chromHMM_enhancer = ifelse(chromHMM == "enhancer", TRUE, FALSE))

setwd(output_dir)

source_barplot <- total_sources %>% 
  ggplot() +
  geom_bar(aes(x = day, group = GENCODE, fill = GENCODE), color = NA, stat = "count", position = "stack") +
  geom_text(stat='count', aes(label=..count.., x = day, y = ..count.. + 4000), size = 6/2.8, color = "black") +
  labs(y = "Number of Transcripts") +
  ylim(0, NA) +
  scale_fill_viridis(discrete = TRUE, option = "C", begin = 0.4, end = 0.7, 
                      labels = c("Novel", "GENCODE"))


pdf("Fig3B_TSS_number.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(source_barplot, width = unit(2.5, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))
dev.off()


upset_df <- total_sources  %>%
  select(id = day, GENCODE, CAGE, chromHMM_prom, chromHMM_enhancer) %>%
  mutate_at(vars(-id), funs(ifelse(. == TRUE, 1, 0)))

for (day in c("d0", "d2", "d4")){

plot_df <- upset_df %>% 
  filter(id == day)
  
pdf(file = paste0("FigS2A-C_", day, "_upsetR_trans_sources.pdf"), onefile = FALSE, height = unit(3, 'cm'), width = unit(4, 'cm')
    ,useDingbats = FALSE)

print(upset(plot_df, order.by = "freq"))
dev.off()
}
