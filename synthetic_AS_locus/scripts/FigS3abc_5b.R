rm(list=ls(all=TRUE))
library(flowCore)
library(openCyto)
library(tidyverse)
library(egg)
library(gridExtra)
library(scales)
library(EnvStats)
library(ggcyto)
library(RColorBrewer)
filter <- dplyr::filter

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


fcs_dir <- "../data/fcs_files_FigS3_2i/"
output_dir <- "../figures/"


flowset <- read.flowSet(path = fcs_dir, min.limit = 0)

setwd(output_dir)

autoplot(flowset, x="FSC-A", y="SSC-A")


rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.35e5, 1.5e5), "SSC-A"= c(0.2e5, 1.5e5))
BoundaryFrame <- Subset (flowset[[1]], rectGate) 

autoplot(flowset[1:12], x="FSC-A", y="SSC-A") + geom_gate(rectGate)


sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = c("FSC-A", "FSC-H"))
singlets_flowset <- Subset (flowset, rectGate %&% sing_g)

autoplot(flowset[1:12], x="FSC-A", y="FSC-H") + geom_gate(sing_g)


extr_frame <- fsApply(singlets_flowset, exprs, simplify = FALSE)
extr_df <- lapply(extr_frame, data.frame)
flow_df <- data.frame(do.call(rbind,extr_df)) %>%
  rownames_to_column("sample") %>%
  transmute(sample, GFP = FITC.A, TagBFP =BV421.A, tdT = PE.A) %>% #sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), 
  separate(sample, c("GenPos", "plasmid", "treatment", "clone"), sep = "_") 

help_df_d2 <- flow_df %>% 
  filter(plasmid == "Empty") 
help_df_gfp <- flow_df %>% 
  filter(plasmid == "Empty") %>% 
  rbind(help_df_d2) %>% 
  mutate(A1 = GFP) %>% 
  select(-GFP, -plasmid, -TagBFP, -tdT, -treatment,-clone) %>% 
  pivot_longer(-1, values_to = "GFP", names_to = "plasmid")
help_df_tagbfp <- flow_df %>% 
  filter(plasmid == "Empty") %>% 
  rbind(help_df_d2) %>% 
  mutate(A1 = TagBFP) %>% 
  select(-GFP, -plasmid, -TagBFP, -tdT, -treatment,-clone) %>% 
  pivot_longer(-1, values_to = "TagBFP", names_to = "plasmid")
help_df_tdt <- flow_df %>% 
  filter(plasmid == "Empty") %>% 
  rbind(help_df_d2) %>% 
  mutate(A1 = tdT) %>% 
  select(-GFP, -plasmid, -TagBFP, -tdT, -treatment,-clone) %>% 
  pivot_longer(-1, values_to = "tdT", names_to = "plasmid")

######### PLOTS

hist_plot <- flow_df %>% 
  dplyr::filter(plasmid != "Empty" & (clone=="C8" | clone=="B6"| clone=="B10")) %>% 
  ggplot() +
  facet_wrap(~clone~treatment, nrow = 6, ncol = 2) +
  geom_density(data = help_df_gfp, 
               aes(x = GFP, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = GFP, ..scaled.., color = treatment)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), limits = c(10, 50000),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#1B9D8E", "#94C3BA"))

fix <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("FigS3a.pdf", fix, width = 21, height = 29.7*2, units = "cm", dpi = 300,useDingbats=FALSE)

hist_plot <- flow_df %>% 
  dplyr::filter(plasmid != "Empty" & (clone=="C8" | clone=="B6"| clone=="B10")) %>% 
  ggplot() +
  facet_wrap(~clone~treatment, nrow = 6, ncol = 2) +
  geom_density(data = help_df_tdt, 
               aes(x = tdT, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = tdT, ..scaled.., color = treatment)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), limits = c(10, 50000),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#EB6D1E", "#F5AC76"))
fix <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("FigS3b.pdf", fix, width = 21, height = 29.7*2, units = "cm", dpi = 300,useDingbats=FALSE)

hist_plot <- flow_df %>% 
  dplyr::filter(plasmid != "Empty" & treatment!="1000nmDox" & (clone=="C8" | clone=="B6"| clone=="B10")) %>% 
  ggplot() +
  facet_wrap(~clone~treatment, nrow = 6, ncol = 1) +
  geom_density(data = help_df_tagbfp, 
               aes(x = TagBFP, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = TagBFP, ..scaled.., color = treatment)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), limits = c(10, 50000),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#473788"))

fix <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("FigS3c.pdf", fix, width = 21, height = 29.7*2, units = "cm", dpi = 300,useDingbats=FALSE)

hist_plot <- flow_df %>% 
  dplyr::filter(plasmid != "Empty" & clone=="C6") %>% 
  ggplot() +
  facet_wrap(~treatment) +
  geom_density(data = help_df_gfp, 
               aes(x = GFP, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = GFP, ..scaled.., color = treatment)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), limits = c(10, 50000),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#1B9D8E", "#94C3BA"))
fix1 <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
hist_plot <- flow_df %>% 
  dplyr::filter(plasmid != "Empty" & clone=="C6") %>% 
  ggplot() +
  facet_wrap(~treatment) +
  geom_density(data = help_df_tdt, 
               aes(x = tdT, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = tdT, ..scaled.., color = treatment)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), limits = c(10, 50000),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#EB6D1E", "#F5AC76"))
fix2 <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
fix<-grid.arrange(fix1,fix2)
ggsave("Fig5b.pdf", fix, dpi = 300,useDingbats=FALSE)
