rm(list=ls(all=TRUE))
library(knitr)
library(flowCore)
library(openCyto)
library(ggcyto)
library(tidyverse)
library(egg)
library(gridExtra)
library(scales)
library(EnvStats)
library(RColorBrewer)
library(pracma)
library("DescTools")
library("ggpubr")

filter <- dplyr::filter

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))


fcs_dir <- "../data/fcs_files_FigS3_differentiation/"
output_dir <- "../figures/"

flowset <- read.flowSet(path = fcs_dir, min.limit = 0)

setwd(output_dir)

autoplot(flowset[1:12], x="FSC-A", y="SSC-A")

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
  transmute(sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), GFP = FITC.A, tdT = PE.A, TBFP = BV421.A) %>%
  separate(sample, c("clone","cellstate", "timepoint", "inicond", "c1","replicate"), sep = "_") 

# duplicate data for 0h timepoint for all c1 concentrations
c1_num <- nrow(flow_df[flow_df$timepoint=="0h",])
c1_vec <- rep(c("0ngDox","2000ngDox"), c1_num)

help_df_0h <- flow_df %>% 
  filter(timepoint=="0h") %>% 
  mutate(rn = row_number()) %>%
  rowwise() %>%
  slice(rep(1, 2)) %>% 
  cbind(c1_new = c1_vec) %>% 
  mutate(c1=c1_new) %>% 
  select(-c1_new,-rn)

flow_df2 <- flow_df %>% 
  filter(timepoint!="0h") %>% 
  rbind(help_df_0h)
rm(help_df_0h)
rm(flow_df)  


# Dataframe for negative data

# duplicate neg data for both inicond
c1_c1<- nrow(flow_df2[flow_df2$clone=="Empty",])
c1_vec <- rep(c("0ngDox","2000ngDox"), c1_c1)
help_df_ini <- flow_df2 %>% 
  filter(clone=="Empty") %>% 
  mutate(rn = row_number()) %>%
  rowwise() %>%
  slice(rep(1, 2)) %>% 
  cbind(c1_new = c1_vec) %>% 
  mutate(c1=c1_new) %>% 
  select(-c1_new,-rn)
# duplicate neg data for all clones
c1_cl <- nrow(help_df_ini[help_df_ini$clone=="Empty",])
c1_vec <- rep(c("SP505-B6","SP505-B10","SP505-C8"), c1_cl)
help_df <- help_df_ini %>% 
  filter(clone=="Empty") %>% 
  mutate(rn = row_number()) %>%
  rowwise() %>%
  slice(rep(1, 3)) %>% 
  cbind(clone_new = c1_vec) %>% 
  mutate(clone=clone_new) %>% 
  select(-clone_new,-rn)
rm(help_df_ini)

stats <- flow_df2 %>% 
  filter(clone != "Empty") %>%
  group_by(clone,cellstate, timepoint, inicond,c1) %>% 
  dplyr::summarize(mfi_gfp = mean(GFP),mfi_tdt = mean(tdT), mfi_tbfp = mean(TBFP)) %>% 
  mutate(timepoint=DigitSum(as.numeric(str_remove_all(timepoint,"[-dh]")))) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 6, 24)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 12, 48)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 9, 72)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 15, 96)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 3, 120)) 
stats_neg <- help_df %>% 
  group_by(clone,cellstate, timepoint, inicond,c1) %>% 
  dplyr::summarize(mfi_gfp_neg = mean(GFP),mfi_tdt_neg = mean(tdT), mfi_tbfp_neg = mean(TBFP)) %>% 
  mutate(timepoint=DigitSum(as.numeric(str_remove_all(timepoint,"[-dh]")))) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 6, 24)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 12, 48)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 9, 72))%>% 
  mutate(timepoint= replace(timepoint, timepoint == 15, 96)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 3, 120)) 
stats_total <- left_join(stats,stats_neg) %>% 
  mutate(mfi_gfp=mfi_gfp-mfi_gfp_neg) %>% 
  mutate(mfi_tdt=mfi_tdt-mfi_tdt_neg) %>%
  mutate(mfi_tbfp=mfi_tbfp-mfi_tbfp_neg) %>%
  select(-mfi_gfp_neg,-mfi_tdt_neg,-mfi_tbfp_neg)

######### PLOTS
mfi_plot_gfp<-stats_total%>% 
  dplyr::filter(clone != "Empty") %>%
  ggplot()+
  geom_point(aes(x=timepoint, y=log10(mfi_gfp),color=clone))+
  geom_line(aes(x=timepoint, y=log10(mfi_gfp),color=clone))+
  facet_wrap(~factor(c1,levels=c("0ngDox","200ngDox","400ngDox","2000ngDox")),ncol=2)
fix1 <- set_panel_size(mfi_plot_gfp, height = unit(2, "cm"), width = unit(1.5, "cm"))
fix<-grid.arrange(fix1)
ggsave(path = output_dir, "FigS3d.pdf", fix, dpi = 300,useDingbats=FALSE)
