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
library(crone)

filter <- dplyr::filter

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

fcs_dir_r1 <- "../data/fcs_files_Fig5/R1/"

fcs_dir_r2 <- "../data/fcs_files_Fig5/R2/"

fcs_dir_r3 <- "../data/fcs_files_Fig5/R3/"

fcs_dir_r4 <- "../data/fcs_files_Fig5/R4/"

fcs_dir_r5 <- "../data/fcs_files_Fig5/R5/"

output_dir <- "../figures/"
setwd(output_dir)

flowset_r1 <- read.flowSet(path = fcs_dir_r1, min.limit = 0)
flowset_r2 <- read.flowSet(path = fcs_dir_r2, min.limit = 0)
flowset_r3 <- read.flowSet(path = fcs_dir_r3, min.limit = 0)
flowset_r4 <- read.flowSet(path = fcs_dir_r4, min.limit = 0)
flowset_r5 <- read.flowSet(path = fcs_dir_r5, min.limit = 0)

autoplot(flowset_r1[1:12], x="FSC-A", y="SSC-A")
autoplot(flowset_r2[1:12], x="FSC-A", y="SSC-A")
autoplot(flowset_r3[1:12], x="FSC-A", y="SSC-A")
autoplot(flowset_r4[1:12], x="FSC-A", y="SSC-A")
autoplot(flowset_r5[1:12], x="FSC-A", y="SSC-A")

rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.35e5, 1.5e5), "SSC-A"= c(0.2e5, 1.5e5))
BoundaryFrame_r1 <- Subset (flowset_r1[[1]], rectGate)
BoundaryFrame_r2 <- Subset (flowset_r2[[1]], rectGate)
BoundaryFrame_r3 <- Subset (flowset_r3[[1]], rectGate) 
BoundaryFrame_r4 <- Subset (flowset_r4[[1]], rectGate) 
BoundaryFrame_r5 <- Subset (flowset_r5[[1]], rectGate) 

autoplot(flowset_r1[1:12], x="FSC-A", y="SSC-A") + geom_gate(rectGate)
autoplot(flowset_r2[1:12], x="FSC-A", y="SSC-A") + geom_gate(rectGate)
autoplot(flowset_r3[1:12], x="FSC-A", y="SSC-A") + geom_gate(rectGate)
autoplot(flowset_r4[1:12], x="FSC-A", y="SSC-A") + geom_gate(rectGate)
autoplot(flowset_r5[1:12], x="FSC-A", y="SSC-A") + geom_gate(rectGate)

sing_g_r1 <- openCyto:::.singletGate(BoundaryFrame_r1, channels = c("FSC-A", "FSC-H"))
sing_g_r2 <- openCyto:::.singletGate(BoundaryFrame_r2, channels = c("FSC-A", "FSC-H"))
sing_g_r3 <- openCyto:::.singletGate(BoundaryFrame_r3, channels = c("FSC-A", "FSC-H"))
sing_g_r4 <- openCyto:::.singletGate(BoundaryFrame_r4, channels = c("FSC-A", "FSC-H"))
sing_g_r5 <- openCyto:::.singletGate(BoundaryFrame_r5, channels = c("FSC-A", "FSC-H"))

singlets_flowset_r1 <- Subset (flowset_r1, rectGate %&% sing_g_r1)
singlets_flowset_r2 <- Subset (flowset_r2, rectGate %&% sing_g_r2)
singlets_flowset_r3 <- Subset (flowset_r3, rectGate %&% sing_g_r3)
singlets_flowset_r4 <- Subset (flowset_r4, rectGate %&% sing_g_r4)
singlets_flowset_r5 <- Subset (flowset_r5, rectGate %&% sing_g_r5)

autoplot(flowset_r1[1:12], x="FSC-A", y="FSC-H") + geom_gate(sing_g_r1)
autoplot(flowset_r2[1:12], x="FSC-A", y="FSC-H") + geom_gate(sing_g_r2)
autoplot(flowset_r3[1:12], x="FSC-A", y="FSC-H") + geom_gate(sing_g_r3)
autoplot(flowset_r4[1:12], x="FSC-A", y="FSC-H") + geom_gate(sing_g_r4)
autoplot(flowset_r5[1:12], x="FSC-A", y="FSC-H") + geom_gate(sing_g_r5)

extr_frame_r1 <- fsApply(singlets_flowset_r1, exprs, simplify = FALSE)
extr_df_r1 <- lapply(extr_frame_r1, data.frame)
flow_df_r1 <- data.frame(do.call(rbind,extr_df_r1)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), GFP = FITC.A, tdT = PE.A, TBFP = BV421.A) %>%
  separate(sample, c("cellstate", "timepoint", "inicond", "c1","replicate"), sep = "_")  %>% 
  mutate(plasmid="SP505")

extr_frame_r2 <- fsApply(singlets_flowset_r2, exprs, simplify = FALSE)
extr_df_r2 <- lapply(extr_frame_r2, data.frame)
flow_df_r2 <- data.frame(do.call(rbind,extr_df_r2)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), GFP = FITC.A, tdT = PE.A, TBFP = BV421.A) %>%
  separate(sample, c("cellstate", "timepoint", "inicond", "c1","replicate"), sep = "_") %>% 
  mutate(plasmid="SP505")

extr_frame_r3 <- fsApply(singlets_flowset_r3, exprs, simplify = FALSE)
extr_df_r3 <- lapply(extr_frame_r3, data.frame)
flow_df_r3 <- data.frame(do.call(rbind,extr_df_r3)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), GFP = FITC.A, tdT = PE.A, TBFP = BV421.A) %>%
  separate(sample, c("cellstate", "timepoint", "inicond", "c1","replicate"), sep = "_") %>% 
  mutate(plasmid="SP505")

extr_frame_r4 <- fsApply(singlets_flowset_r4, exprs, simplify = FALSE)
extr_df_r4 <- lapply(extr_frame_r4, data.frame)
flow_df_r4 <- data.frame(do.call(rbind,extr_df_r4)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), GFP = FITC.A, tdT = PE.A, TBFP = BV421.A) %>%
  separate(sample, c("cellstate", "timepoint", "inicond", "c1","replicate"), sep = "_") %>% 
  mutate(plasmid="SP505")

extr_frame_r5 <- fsApply(singlets_flowset_r5, exprs, simplify = FALSE)
extr_df_r5 <- lapply(extr_frame_r5, data.frame)
flow_df_r5 <- data.frame(do.call(rbind,extr_df_r5)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub("_[A-H][0-9].*.fcs.*", "", sample), GFP = FITC.A, tdT = PE.A, TBFP = BV421.A) %>%
  separate(sample, c("plasmid", "inicond", "timepoint", "c1","replicate"), sep = "_") %>% 
  mutate(plasmid= replace(plasmid, plasmid == "A4", "SP538")) %>% 
  mutate(plasmid= replace(plasmid, plasmid == "C6", "SP505")) %>%
  mutate(inicond= replace(inicond, inicond == "ini1", "Ini1-noDox")) %>%
  mutate(inicond= replace(inicond, inicond == "ini2", "Ini2-highDox")) %>%
  mutate(replicate= replace(replicate, replicate == "R1", "R5"), replicate= replace(replicate, replicate == "R2", "R6"), replicate= replace(replicate, replicate == "R3", "R7")) %>%
  mutate(c1= replace(c1, c1 == "0", "0ngDox"),c1= replace(c1, c1 == "200", "200ngDox"),c1= replace(c1, c1 == "400", "400ngDox"),c1= replace(c1, c1 == "2000", "2000ngDox")) %>% 
  mutate(c1= replace(c1, c1 == " ", "noc1"),c1= replace(c1, c1 == "  ", "noc1")) %>%
  mutate(inicond= replace(inicond, plasmid == "LP", "Neg")) %>%
  mutate(cellstate="2i") %>% 
  mutate_at(c('timepoint'), as.numeric)

# Combine replicates and filter the data 
flow_df <- rbind(flow_df_r1,flow_df_r2,flow_df_r3,flow_df_r4, flow_df_r5) %>% 
  filter(cellstate!="EpiSC") %>% 
  filter(c1!="10ngDox"&c1!="50ngDox"&c1!="100ngDox"&c1!="150ngDox"&c1!="250ngDox"&c1!="300ngDox"&c1!="800ngDox"&c1!="1000ngDox") %>% 
  mutate(cellstate= replace(cellstate, cellstate == "ESC", "2i")) %>% 
  mutate(replicate= replace(replicate, replicate == "R1", "R2")) %>% 
  mutate(replicate= replace(replicate, replicate == "R2", "R1")) %>% 
  mutate(replicate= replace(replicate, replicate == "R3", "R2")) %>%
  mutate(replicate= replace(replicate, replicate == "R4", "R3")) %>% 
  mutate(replicate= replace(replicate, replicate == "R5", "R4")) %>% 
  mutate(replicate= replace(replicate, replicate == "R6", "R5")) %>% 
  mutate(replicate= replace(replicate, replicate == "R7", "R6"))
rm(flow_df_r1,flow_df_r2,flow_df_r3,flow_df_r4,flow_df_r5)

# duplicate data for 0h timepoint for all c1 concentrations
c1_num <- nrow(flow_df[flow_df$timepoint=="0h"|flow_df$timepoint=="0",])
c1_vec <- rep(c("0ngDox","200ngDox",
                "400ngDox","2000ngDox"), c1_num)

help_df_0h <- flow_df %>% 
  filter(timepoint=="0h"|timepoint=="0") %>% 
  mutate(rn = row_number()) %>%
  rowwise() %>%
  slice(rep(1, 8)) %>% 
  cbind(c1_new = c1_vec) %>% 
  mutate(c1=c1_new) %>% 
  select(-c1_new,-rn)

flow_df2 <- flow_df %>% 
  filter(c1!="noc1") %>% 
  rbind(help_df_0h)
rm(help_df_0h)
rm(flow_df)

# Dataframe for negative data

help_df1_1 <- flow_df2 %>% 
  filter(inicond == "Neg" & plasmid=="LP") %>% 
  mutate(inicond = "Ini1-noDox") %>% 
  mutate(plasmid="SP538")
help_df1_2 <- flow_df2 %>% 
  filter(inicond == "Neg") %>% 
  mutate(inicond = "Ini1-noDox") %>% 
  mutate(plasmid= replace(plasmid, plasmid == "LP", "SP505"))
help_df2_1 <- flow_df2 %>% 
  filter(inicond == "Neg" & plasmid=="LP") %>% 
  mutate(inicond = "Ini2-highDox") %>% 
  mutate(plasmid="SP538")
help_df2_2 <- flow_df2 %>% 
  filter(inicond == "Neg") %>% 
  mutate(inicond = "Ini2-highDox") %>% 
  mutate(plasmid= replace(plasmid, plasmid == "LP", "SP505"))
help_df_ic <- rbind(help_df1_1,help_df1_2,help_df2_1,help_df2_2)
help_df_c1_1 <- help_df_ic %>% 
  mutate(c1 ="0ngDox")
help_df_c1_2 <- help_df_ic %>% 
  mutate(c1 ="200ngDox")
help_df_c1_3 <- help_df_ic %>% 
  mutate(c1 ="400ngDox")
help_df_c1_4 <- help_df_ic %>% 
  mutate(c1 ="2000ngDox")

help_df <- rbind(help_df_c1_1,help_df_c1_2,help_df_c1_3,help_df_c1_4) 
help_df_2ilw <- help_df %>% 
  filter(cellstate=="2i-LIFw")
help_df_lw <- help_df %>% 
  filter(cellstate=="LIFw")
rm(help_df_c1_1,help_df_c1_2,help_df_c1_3,help_df_c1_4)

# Print minimal number of live single cell ESC events acquired per sample
events<-flow_df2 %>% 
  filter(inicond != "Neg") %>%
  group_by(cellstate, timepoint, inicond,c1,replicate,plasmid) %>%
  count(cellstate, timepoint, inicond,c1,replicate,plasmid, sort=TRUE)
tail(events,n=10)

stats <- flow_df2 %>% 
  filter(inicond != "Neg") %>%
  filter(replicate != "R6") %>% #too few cells
  group_by(cellstate, timepoint, inicond,c1,replicate,plasmid) %>% 
  dplyr::summarize(mfi_gfp = mean(GFP),mfi_tdt = mean(tdT), mfi_tbfp = mean(TBFP)) %>% 
  mutate(timepoint=DigitSum(as.numeric(str_remove_all(timepoint,"[-dh]")))) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 6, 24)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 12, 48)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 9, 72)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 15, 96)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 3, 120)) %>% 
  mutate(replicate= replace(replicate, replicate == "R4" & plasmid=="SP538", "R1")) %>% 
  mutate(replicate= replace(replicate, replicate == "R5" & plasmid=="SP538", "R2"))  
stats_neg <- help_df %>% 
  filter(replicate != "R6") %>% #too few cells
  group_by(cellstate, timepoint, inicond,c1,replicate,plasmid) %>% 
  dplyr::summarize(mfi_gfp_neg = mean(GFP),mfi_tdt_neg = mean(tdT), mfi_tbfp_neg = mean(TBFP)) %>% 
  mutate(timepoint=DigitSum(as.numeric(str_remove_all(timepoint,"[-dh]")))) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 6, 24)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 12, 48)) %>%
  mutate(timepoint= replace(timepoint, timepoint == 9, 72)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 15, 96)) %>% 
  mutate(timepoint= replace(timepoint, timepoint == 3, 120)) %>% 
  mutate(replicate= replace(replicate, replicate == "R4" & plasmid=="SP538", "R1")) %>% 
  mutate(replicate= replace(replicate, replicate == "R5" & plasmid=="SP538", "R2")) 
stats_total <- left_join(stats,stats_neg) %>% 
  mutate(mfi_gfp=mfi_gfp-mfi_gfp_neg) %>% 
  mutate(mfi_tdt=mfi_tdt-mfi_tdt_neg) %>%
  mutate(mfi_tbfp=mfi_tbfp-mfi_tbfp_neg) %>%
  select(-mfi_gfp_neg,-mfi_tdt_neg,-mfi_tbfp_neg)
stats_total[stats_total < 0] <- 0 # set negative values to 0
# Create a group-means data set,
mean_reps <-stats_total %>% 
  filter(timepoint<=96) %>% 
  group_by(cellstate,timepoint,inicond,c1, plasmid) %>% 
  summarise(mfi_gfp_rep = mean(mfi_gfp),mfi_tdt_rep = mean(mfi_tdt), mfi_tbfp_rep = mean(mfi_tbfp)) 

# create dataframe containing difference in MFI between 2 initial conditions
diff_inis <- stats_total %>% 
  dplyr::filter(inicond != "Neg") %>%
  select(-mfi_tdt,-mfi_tbfp) %>% 
  pivot_wider(names_from=inicond,values_from=mfi_gfp) %>% 
  mutate(mfi_gfp_diff = `Ini1-noDox`-`Ini2-highDox`) %>% 
  mutate(mfi_gfp_diff_log = log10(`Ini1-noDox`)-log10(`Ini2-highDox`)) %>% 
  select(-`Ini1-noDox`,-`Ini2-highDox`)
mean_diff_inis <- diff_inis %>% 
  filter(timepoint<=96) %>% 
  group_by(cellstate,timepoint,c1,plasmid) %>% 
  summarise(mfi_gfp_diff_rep = mean(mfi_gfp_diff),mfi_gfp_diff_log_rep = mean(mfi_gfp_diff_log))

diff_inis_tdt <- stats_total %>% 
  dplyr::filter(inicond != "Neg") %>%
  select(-mfi_gfp,-mfi_tbfp) %>% 
  pivot_wider(names_from=inicond,values_from=mfi_tdt) %>% 
  mutate(mfi_tdt_diff = `Ini1-noDox`-`Ini2-highDox`) %>% 
  mutate(mfi_tdt_diff_log = log10(`Ini1-noDox`)-log10(`Ini2-highDox`)) %>% 
  select(-`Ini1-noDox`,-`Ini2-highDox`)
mean_diff_inis_tdt <- diff_inis_tdt %>% 
  filter(timepoint<=96) %>% 
  filter(!mfi_tdt_diff_log==-Inf) %>% #remove values with bdg corrected FI = 0 from calculation of mean to prevent -Inf
  group_by(cellstate,timepoint,c1,plasmid) %>% 
  summarise(mfi_tdt_diff_rep = mean(mfi_tdt_diff),mfi_tdt_diff_log_rep = mean(mfi_tdt_diff_log))

# Calculate FC to each respective non-Dox treated timepoint (0ng Dox and coming from Ini1)
c1_num <- nrow(stats_total[(stats_total$c1=="0ngDox"&stats_total$inicond=="Ini1-noDox"),])
c1_vec <- rep(c("0ngDox","200ngDox",
                "400ngDox","2000ngDox"), c1_num)

help_stats_nodox<- stats_total %>% 
  filter(c1=="0ngDox"&inicond=="Ini1-noDox") %>% 
  mutate(rn = row_number()) %>%
  rowwise() %>%
  slice(rep(1, 8)) %>% 
  cbind(c1_new = c1_vec) %>% 
  mutate(c1=c1_new, mfi_gfp_nodox=mfi_gfp, mfi_tdt_nodox=mfi_tdt, mfi_tbfp_nodox=mfi_tbfp) %>% 
  select(-c1_new,-rn, -mfi_gfp,-mfi_tdt,-mfi_tbfp)

ini_num <- nrow(help_stats_nodox[help_stats_nodox$inicond=="Ini1-noDox",])
ini_vec <- rep(c("Ini1-noDox","Ini2-highDox"), ini_num)

help_stats_nodox2 <- help_stats_nodox %>% 
  mutate(rn = row_number()) %>%
  rowwise() %>%
  slice(rep(1, 2)) %>% 
  cbind(ini_new = ini_vec) %>% 
  mutate(inicond=ini_new) %>% 
  select(-ini_new)

fc_total <- left_join(stats_total,help_stats_nodox2) %>% 
  mutate(fc_gfp=mfi_gfp/mfi_gfp_nodox) %>% 
  mutate(fc_tdt=mfi_tdt/mfi_tdt_nodox) %>%
  mutate(fc_tbfp=mfi_tbfp/mfi_tbfp_nodox) %>%
  select(-mfi_gfp_nodox,-mfi_tdt_nodox,-mfi_tbfp_nodox)

mean_fc_reps <-fc_total %>% 
  filter(timepoint<=96) %>% 
  group_by(cellstate,timepoint,inicond,c1,plasmid) %>% 
  summarise(mfi_gfp_rep = mean(mfi_gfp),mfi_tdt_rep = mean(mfi_tdt), mfi_tbfp_rep = mean(mfi_tbfp), mfc_gfp_rep=mean(fc_gfp), mfc_tdt_rep=mean(fc_tdt),
            mfc_tbfp_rep=mean(fc_tbfp))

### Calculate ratio of MFI in Ini1 to Ini2
help_stats_ini1 <- stats_total %>% 
  filter(inicond=="Ini1-noDox") %>% 
  group_by(cellstate, timepoint,c1,replicate,plasmid) %>% 
  mutate(gfp_ini1=mfi_gfp,tdt_ini1=mfi_tdt,tbfp_ini1=mfi_tbfp) %>% 
  select(-mfi_gfp,-mfi_tdt,-mfi_tbfp,-inicond)
help_stats_ini2 <- stats_total %>% 
  filter(inicond=="Ini2-highDox") %>% 
  group_by(cellstate, timepoint,c1,replicate,plasmid) %>% 
  mutate(gfp_ini2=mfi_gfp,tdt_ini2=mfi_tdt,tbfp_ini2=mfi_tbfp) %>% 
  select(-mfi_gfp,-mfi_tdt,-mfi_tbfp,-inicond)

mfi_rat <- left_join(help_stats_ini1,help_stats_ini2) %>% 
  mutate(rat_gfp=gfp_ini1/gfp_ini2,rat_tdt=tdt_ini1/tdt_ini2,rat_tbfp=tbfp_ini1/tbfp_ini2) 

###### SIMULATION of GFP and tdT kinetics in the absence of memory
# parameters
kdeg_gfp <- log(2)/(26*0.59) # literature half life GFP:26h, fold change of 0.59 when comparing (GFP-FKBP-1uM-Shld)/(WT-GFP) Expression levels 
kdeg_tdt <- log(2)/(26*0.86)
kdil <- log(2)/12 # dilution due to cell division (cell cycle length ~11h (Waisman2018))
kdd_gfp <- kdeg_gfp + kdil # total rate of GFP removal due to degradation & dilution
kdd_tdt <- kdeg_tdt + kdil # total rate of GFP removal due to degradation & dilution
time <- (0:100)

### Analytical solution of ODE 2i
gfp_ode1 <- (1-exp(time*(-kdd_gfp)))# high Dox -> no Dox with GFP(t=0)=0 and kprod=kdeg
gfp_ode2 <- (exp(time*(-kdd_gfp)))# no Dox -> high Dox with GFP(t=0)=1 and kprod=0
tdt_ode1 <- (exp(time*(-kdd_tdt)))# high Dox -> no Dox
tdt_ode2 <- (1-exp(time*(-kdd_tdt)))# no Dox -> high Dox


### Scaling of simulation to FACS data 
mean_sc <- mean_reps %>% 
  filter(((inicond=="Ini1-noDox" & c1=="0ngDox")|(inicond=="Ini2-highDox" & c1=="2000ngDox"))&cellstate!="2i-LIFw") %>% 
  group_by(timepoint,cellstate,inicond,plasmid) %>% 
  dplyr::summarize(gfp = mean(mfi_gfp_rep),tdt = mean(mfi_tdt_rep)) %>%
  pivot_wider(names_from=c(cellstate,inicond,plasmid),values_from = c(gfp,tdt))

# 2i: scale to time-average
gfp_ode1_sc_2i_505 <- gfp_ode1*(mean(mean_sc$'gfp_2i_Ini1-noDox_SP505')-mean(mean_sc$'gfp_2i_Ini2-highDox_SP505'))+mean(mean_sc$'gfp_2i_Ini2-highDox_SP505')
gfp_ode2_sc_2i_505 <- gfp_ode2*(mean(mean_sc$'gfp_2i_Ini1-noDox_SP505')-mean(mean_sc$'gfp_2i_Ini2-highDox_SP505'))+mean(mean_sc$'gfp_2i_Ini2-highDox_SP505')
gfp_ode1_sc_2i_538 <- gfp_ode1*(mean(mean_sc$'gfp_2i_Ini1-noDox_SP538',na.rm=TRUE)-mean(mean_sc$'gfp_2i_Ini2-highDox_SP538',na.rm=TRUE))+mean(mean_sc$'gfp_2i_Ini2-highDox_SP538',na.rm=TRUE)
gfp_ode2_sc_2i_538 <- gfp_ode2*(mean(mean_sc$'gfp_2i_Ini1-noDox_SP538',na.rm=TRUE)-mean(mean_sc$'gfp_2i_Ini2-highDox_SP538',na.rm=TRUE))+mean(mean_sc$'gfp_2i_Ini2-highDox_SP538',na.rm=TRUE)

tdt_ode1_sc_2i_505 <- tdt_ode1*(mean(mean_sc$'tdt_2i_Ini2-highDox_SP505')-mean(mean_sc$'tdt_2i_Ini1-noDox_SP505'))+mean(mean_sc$'tdt_2i_Ini1-noDox_SP505')
tdt_ode2_sc_2i_505 <- tdt_ode2*(mean(mean_sc$'tdt_2i_Ini2-highDox_SP505')-mean(mean_sc$'tdt_2i_Ini1-noDox_SP505'))+mean(mean_sc$'tdt_2i_Ini1-noDox_SP505')
tdt_ode1_sc_2i_538 <- tdt_ode1*(mean(mean_sc$'tdt_2i_Ini2-highDox_SP538',na.rm=TRUE)-mean(mean_sc$'tdt_2i_Ini1-noDox_SP538',na.rm=TRUE))+mean(mean_sc$'tdt_2i_Ini1-noDox_SP538',na.rm=TRUE)
tdt_ode2_sc_2i_538 <- tdt_ode2*(mean(mean_sc$'tdt_2i_Ini2-highDox_SP538',na.rm=TRUE)-mean(mean_sc$'tdt_2i_Ini1-noDox_SP538',na.rm=TRUE))+mean(mean_sc$'tdt_2i_Ini1-noDox_SP538',na.rm=TRUE)

gfp_ode<- bind_cols(time,gfp_ode1_sc_2i_505,gfp_ode2_sc_2i_505,gfp_ode1_sc_2i_538,gfp_ode2_sc_2i_538) %>% 
  transmute(timepoint=...1, gfp_ode1_2i_505=...2, gfp_ode2_2i_505=...3,gfp_ode1_2i_538=...4, gfp_ode2_2i_538=...5) %>% 
  pivot_longer(!timepoint,names_to = c("reporter","inicond","cellstate","plasmid") , values_to = c("mfi_gfp_ode"), names_sep="_") %>% 
  mutate(c1="0ngDox", c1=replace(c1,inicond == "ode2","2000ngDox")) %>% 
  mutate(inicond = ifelse(inicond == 'ode2', 'Ini1-noDox', 'Ini2-highDox')) %>% 
  mutate(plasmid = ifelse(plasmid == '505', 'SP505', 'SP538')) %>% 
  select(-reporter)

tdt_ode<- bind_cols(time,tdt_ode1_sc_2i_505,tdt_ode2_sc_2i_505,tdt_ode1_sc_2i_538,tdt_ode2_sc_2i_538) %>% 
  transmute(timepoint=...1, tdt_ode1_2i_505=...2, tdt_ode2_2i_505=...3,tdt_ode1_2i_538=...4, tdt_ode2_2i_538=...5) %>% 
  pivot_longer(!timepoint,names_to = c("reporter","inicond","cellstate","plasmid") , values_to = c("mfi_tdt_ode"), names_sep="_") %>% 
  mutate(c1="0ngDox", c1=replace(c1,inicond == "ode2","2000ngDox")) %>% 
  mutate(inicond = ifelse(inicond == 'ode2', 'Ini1-noDox', 'Ini2-highDox')) %>% 
  mutate(plasmid = ifelse(plasmid == '505', 'SP505', 'SP538')) %>% 
  select(-reporter)

##### Quantify memory by testing different offsets for the ode simulation:
# With which offset does ODE simulation best explain the data?

# define range of ode offsets (=memory of the system)
offset<-cbind(0:48)
# measurement timepoints
tp<- unique(stats_total$timepoint)
tp<- tp[tp <= 96]
tp<-sort(tp)

ssqr_gfp = data.frame()#dataframe to save sum of squared residuals for each offset
ssqr_tdt = data.frame()
# loop over offsets
for (o in offset) {
  # shift ode simulation by offset
  gfp_ode_shift <- gfp_ode %>% 
    mutate(timepoint = timepoint + offset(o)) 
  tdt_ode_shift <- tdt_ode %>% 
    mutate(timepoint = timepoint + offset(o)) 
  # select rows from ode simulation that correspond to measurement timepoints
  gfp_ode_shift <- gfp_ode_shift %>% filter(timepoint %in% tp)
  tdt_ode_shift <- tdt_ode_shift %>% filter(timepoint %in% tp)
  # add measurement timpoints that fall into window 0:offset-1 with initial gfp value
  tp_add <- tp[tp < o]
  for (j in tp_add) {
    gfp_ode_shift<- gfp_ode_shift %>% 
    bind_rows(gfp_ode %>% filter(timepoint == 0) %>% mutate(timepoint =j))
    tdt_ode_shift<- tdt_ode_shift %>% 
      bind_rows(tdt_ode %>% filter(timepoint == 0) %>% mutate(timepoint =j))
  }
  # calculate SQRs between shifted ODE and data
  data_gfp<-stats_total %>% 
    filter((inicond=="Ini2-highDox"&c1=="0ngDox")|(inicond=="Ini1-noDox"&c1=="2000ngDox")) %>% 
    group_by(cellstate,timepoint,inicond,c1,plasmid,replicate) %>% 
    summarize(mean_mfi_gfp = mean(mfi_gfp))
  sqr_gfp<-left_join(gfp_ode_shift, data_gfp, by = c("timepoint", "inicond","plasmid","cellstate", "c1")) %>% 
    mutate(sqr=(mfi_gfp_ode-mean_mfi_gfp)^2) %>% 
    mutate(offset=o) %>% 
    filter(!is.na(mean_mfi_gfp)) %>% 
    group_by(cellstate,inicond,c1,plasmid,replicate) %>% 
    summarize(sum_sqr=sum(sqr),offset=mean(offset))
  ssqr_gfp <- rbind(ssqr_gfp,sqr_gfp)
  
  data_tdt<-stats_total %>% 
    filter((inicond=="Ini2-highDox"&c1=="0ngDox")|(inicond=="Ini1-noDox"&c1=="2000ngDox")) %>% 
    group_by(cellstate,timepoint,inicond,c1,plasmid,replicate) %>% 
    summarize(mean_mfi_tdt = mean(mfi_tdt))
  sqr_tdt<-left_join(tdt_ode_shift, data_tdt, by = c("timepoint", "inicond","plasmid","cellstate", "c1")) %>% 
    mutate(sqr=(mfi_tdt_ode-mean_mfi_tdt)^2) %>% 
    mutate(offset=o) %>% 
    filter(!is.na(mean_mfi_tdt)) %>% # remove rows with missing measurements (SP538 not measured at 8h)
    group_by(cellstate,inicond,c1,plasmid,replicate) %>% 
    summarize(sum_sqr=sum(sqr),offset=mean(offset))
  ssqr_tdt <- rbind(ssqr_tdt,sqr_tdt)
}
# select minimal ssqr
min_ssqr_gfp <- ssqr_gfp %>% group_by(cellstate, inicond, c1,plasmid, replicate) %>%  
  filter(sum_sqr == min(sum_sqr)) %>% 
  group_by(cellstate, inicond, c1,plasmid) %>%  
  summarise(mean_offset = mean(offset, na.rm = TRUE), sd_offset = sd(offset, na.rm = TRUE), min_offset=min(offset, na.rm = TRUE), max_offset = max(offset, na.rm = TRUE), sum_sqr=mean(sum_sqr, na.rm = TRUE))  
min_ssqr_tdt <- ssqr_tdt %>% group_by(cellstate, inicond, c1,plasmid, replicate) %>%  
  filter(sum_sqr == min(sum_sqr)) %>% 
  group_by(cellstate, inicond, c1,plasmid) %>%  
  summarise(mean_offset = mean(offset, na.rm = TRUE), sd_offset = sd(offset, na.rm = TRUE), , min_offset=min(offset, na.rm = TRUE), max_offset = max(offset, na.rm = TRUE), sum_sqr=mean(sum_sqr, na.rm = TRUE)) 


# to include shifted ode in plot, construct tibble with optimally shifted ode
gfp_ode_memory <- gfp_ode %>% 
  left_join(min_ssqr_gfp,by = c( "inicond","plasmid","cellstate", "c1")) %>% 
  mutate(timepoint = timepoint + mean_offset) 
for (i in unique(gfp_ode_memory$plasmid)) {
  for(j in unique(gfp_ode_memory$inicond)) {
    for(s in unique(gfp_ode_memory$cellstate)) {
add_gfp_ode_memory <- gfp_ode_memory %>%
  filter(timepoint==mean_offset&plasmid==i&inicond==j&cellstate==s) %>%
  group_by(timepoint,inicond,cellstate,plasmid,mfi_gfp_ode,c1,sum_sqr,mean_offset) %>%
  slice(rep(1, mean_offset)) %>%
  ungroup() %>%
  mutate(timepoint=row_number())
  gfp_ode_memory<-rbind(gfp_ode_memory,add_gfp_ode_memory)
    }
  }
}

tdt_ode_memory <- tdt_ode %>% 
  left_join(min_ssqr_tdt,by = c( "inicond","plasmid","cellstate", "c1")) %>% 
  mutate(timepoint = timepoint + mean_offset) 
for (i in unique(tdt_ode_memory$plasmid)) {
  for(j in unique(tdt_ode_memory$inicond)) {
    for(s in unique(tdt_ode_memory$cellstate)) {
      add_tdt_ode_memory <- tdt_ode_memory %>%
        filter(timepoint==mean_offset&plasmid==i&inicond==j&cellstate==s) %>%
        group_by(timepoint,inicond,cellstate,plasmid,mfi_tdt_ode,c1,sum_sqr,mean_offset) %>%
        slice(rep(1, mean_offset)) %>%
        ungroup() %>%
        mutate(timepoint=row_number())
      tdt_ode_memory<-rbind(tdt_ode_memory,add_tdt_ode_memory)
    }
  }
}


########## Test for statistical significance: Population GFP and tdt means different between 2 inis?

stat_gfp<-stats_total %>% 
  filter(timepoint<=96) %>% 
  select(-mfi_tdt,-mfi_tbfp) %>% 
  pivot_wider(names_from = inicond, values_from = mfi_gfp) %>% 
  group_by(cellstate,timepoint,c1,plasmid) %>% 
  summarise(pval.2sample_gfp=t.test(log10(`Ini1-noDox`),log10(`Ini2-highDox`),paired=F)$p.value)

stat_tdt<-stats_total %>% 
  filter(timepoint<=96) %>% 
  mutate(mfi_tdt=replace(mfi_tdt, mfi_tdt == 0, 0.00001)) %>% 
  select(-mfi_gfp,-mfi_tbfp) %>% 
  pivot_wider(names_from = inicond, values_from = mfi_tdt) %>% 
  group_by(cellstate,timepoint,c1,plasmid) %>% 
  summarise(pval.2sample_tdt=t.test(log10(`Ini1-noDox`),log10(`Ini2-highDox`),paired=F)$p.value)

stats_total.stat <- stats_total %>% left_join(stat_gfp) %>% left_join(stat_tdt) 
# Create dataframe to add p-value * to plots
sel.stat = stats_total.stat %>% 
  filter(replicate=="R1") %>% 
  mutate(y_gfp_2i=3.6,y_gfp_lifw=4.1,y_tdt_2i=4.25,y_tdt_lifw=5.15,y_gfp_2i_538=4.23,y_tdt_2i_538=5,lab_gfp=ifelse(signif(pval.2sample_gfp,2)<0.05,'*',''),lab_tdt=ifelse(signif(pval.2sample_tdt,2)<0.05,'*',''))



#### PLOTS
s_m = T # set to true to include simulation with optimal shift (=memory)

##### Figure 5D
# MFI GFP
mfi_gfp_2i_sp505<-stats_total%>% 
  dplyr::filter(inicond != "Neg") %>%
  filter(timepoint<=96 & cellstate=="2i" & plasmid=="SP505") %>% 
  ggplot()+
  geom_point(aes(x=timepoint, y=log10(mfi_gfp),color=inicond, shape=inicond),size = 1, alpha = 2/3)+ 
  geom_line(data=mean_reps %>% filter(timepoint<=96& cellstate=="2i" & plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_gfp_rep),color=inicond))+
  # Add simulated GFP kinetics with optimal shift = memory (least ssqr)
  {if(s_m)geom_line(data=gfp_ode_memory %>% filter(cellstate=="2i"&plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_gfp_ode)),colour="black",size = 0.5,linetype = 1)}+
  # plot statistical test for difference between 2 initial conditions
  geom_text(data=sel.stat%>% filter(timepoint<=96& cellstate=="2i"& plasmid=="SP505"),aes(x=timepoint,y=y_gfp_2i, label=lab_gfp), size=5, color='black')  +
  scale_fill_manual(values=c("#1B9D8E", "#94C3BA"))+
  scale_shape_manual(values = c(16,17,18)) + # Customize point symbols
  scale_color_manual(values=c("#1B9D8E", "#94C3BA"))+
  scale_x_continuous(breaks=c(0,8,24,48,72,96), labels=c(0,0.33,1,2,3,4), expand = c(0, 0), limits=c(-4,100))+
  scale_y_continuous(limits=c(2.5,3.75),breaks=c(2.5, 2.75,3,3.25,3.5,3.75), name='MFI GFP [log10]') +
  facet_wrap(cellstate~factor(c1,levels=c("0ngDox","200ngDox",
                                          "400ngDox", "2000ngDox")),ncol=4)
fix_mfi_gfp_2i_sp505 <- set_panel_size(mfi_gfp_2i_sp505, height = unit(2, "cm"), width = unit(1.5, "cm"))

# MFI tdTomato
mfi_tdt_2i_sp505<-stats_total%>% 
  dplyr::filter(inicond != "Neg") %>%
  filter(timepoint<=96 & cellstate=="2i" & plasmid=="SP505") %>% 
  ggplot()+
  geom_point(aes(x=timepoint, y=log10(mfi_tdt), color=inicond, shape=inicond),size = 1, alpha = 2/3)+ 
  geom_line(data=mean_reps %>% filter(timepoint<=96& cellstate=="2i"& plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_tdt_rep),color=inicond))+
  # Add simulated tdt kinetics
  {if(s_m)geom_line(data=tdt_ode_memory %>% filter(cellstate=="2i"&plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_tdt_ode)),colour="black",size = 0.5,linetype = 1)}+ 
  geom_text(data=sel.stat%>% filter(timepoint<=96& cellstate=="2i"& plasmid=="SP505"),aes(x=timepoint,y=y_tdt_2i, label=lab_tdt), size=5, color='black')  +
  scale_fill_manual(values=c("#EB6D1E", "#F5AC76"))+
  scale_shape_manual(values = c(16,17)) +
  scale_color_manual(values=c("#EB6D1E", "#F5AC76"))+
  scale_x_continuous(breaks=c(0,8,24,48,72,96), labels=c(0,0.33,1,2,3,4), expand = c(0, 0), limits=c(-4,100))+
  scale_y_continuous(limits=c(-1,4.75),breaks=c((2.5)-1,0,1,2,3,4), name='MFI tdT [log10]') +
  facet_wrap(cellstate~factor(c1,levels=c("0ngDox","200ngDox",
                                          "400ngDox", "2000ngDox")),ncol=4)
fix_mfi_tdt_2i_sp505<- set_panel_size(mfi_tdt_2i_sp505, height = unit(2, "cm"), width = unit(1.5, "cm"))

fig5d<-grid.arrange(fix_mfi_gfp_2i_sp505,fix_mfi_tdt_2i_sp505)
ggsave(path = output_dir, "Fig5d.pdf", fig5d, dpi = 300,
       useDingbats=FALSE)

#### Figure 5E
#MFI GFP
mfi_gfp_lifw_sp505<-stats_total%>% 
  dplyr::filter(inicond != "Neg") %>%
  filter(timepoint<=96 & cellstate=="LIFw" & plasmid=="SP505") %>% 
  ggplot()+
  geom_point(aes(x=timepoint, y=log10(mfi_gfp),color=inicond, shape=inicond),size = 1, alpha = 2/3)+ 
  geom_line(data=mean_reps %>% filter(timepoint<=96& cellstate=="LIFw" & plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_gfp_rep),color=inicond))+
  # Add simulated GFP kinetics
  {if(s_m)geom_line(data=gfp_ode_memory %>% filter(cellstate=="LIFw"&plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_gfp_ode)),colour="grey",size = 0.5,linetype = 1)}+ 
   geom_text(data=sel.stat%>% filter(timepoint<=96& cellstate=="LIFw"& plasmid=="SP505"),aes(x=timepoint,y=y_gfp_lifw, label=lab_gfp), size=5, color='black')  +
  scale_fill_manual(values=c("#1B9D8E", "#94C3BA"))+
  scale_shape_manual(values = c(16,17)) +
  scale_color_manual(values=c("#1B9D8E", "#94C3BA"))+ 
  scale_x_continuous(breaks=c(0,8,24,48,72,96), labels=c(0,0.33,1,2,3,4), expand = c(0, 0), limits=c(-4,100))+
  scale_y_continuous(limits=c(1.75,4.35),breaks=c(2,2.5,3,3.5,4), name='MFI GFP [log10]') +
  facet_wrap(cellstate~factor(c1,levels=c("0ngDox","200ngDox",
                                          "400ngDox", "2000ngDox")),ncol=4)
fix_mfi_gfp_lifw_sp505 <- set_panel_size(mfi_gfp_lifw_sp505, height = unit(2, "cm"), width = unit(1.5, "cm"))

#MFI tdT
mfi_tdt_lifw_sp505<-stats_total%>% 
  dplyr::filter(inicond != "Neg") %>%
  filter(timepoint<=96 & cellstate=="LIFw" & plasmid=="SP505") %>% 
  ggplot()+
  geom_point(aes(x=timepoint, y=log10(mfi_tdt),color=inicond, shape=inicond),size = 1, alpha = 2/3)+
  geom_line(data=mean_reps %>% filter(timepoint<=96& cellstate=="LIFw"& plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_tdt_rep),color=inicond))+
  # Add simulated tdt kinetics
  {if(s_m)geom_line(data=tdt_ode_memory %>% filter(cellstate=="LIFw"&plasmid=="SP505"),aes(x=timepoint, y=log10(mfi_tdt_ode)),colour="grey",size = 0.5,linetype = 1)}+ 
  geom_text(data=sel.stat%>% filter(timepoint<=96& cellstate=="LIFw"& plasmid=="SP505"),aes(x=timepoint,y=y_tdt_lifw, label=lab_tdt), size=5, color='black')  +
  geom_text(data=sel.stat%>% filter(timepoint<=96& cellstate=="LIFw"& plasmid=="SP505"),aes(x=timepoint,y=y_tdt_lifw, label=lab_tdt), size=5, color='black')  +
  scale_fill_manual(values=c("#EB6D1E", "#F5AC76"))+
  scale_shape_manual(values = c(16,17)) +
  scale_color_manual(values=c("#EB6D1E", "#F5AC76"))+
  scale_x_continuous(breaks=c(0,8,24,48,72,96), labels=c(0,0.33,1,2,3,4), expand = c(0, 0), limits=c(-4,100))+
  scale_y_continuous(limits=c(0.5,5.5),breaks=c(1,2,3,4,5), name='MFI tdT [log10]') +
  facet_wrap(cellstate~factor(c1,levels=c("0ngDox","200ngDox",
                                          "400ngDox", "2000ngDox")),ncol=4)
fix_mfi_tdt_lifw_sp505 <- set_panel_size(mfi_tdt_lifw_sp505, height = unit(2, "cm"), width = unit(1.5, "cm"))

fig5e<-grid.arrange(fix_mfi_gfp_lifw_sp505,fix_mfi_tdt_lifw_sp505)
ggsave(path = output_dir, "Fig5e.pdf", fig5e, dpi = 300,
       useDingbats=FALSE)


# Delta MFI GFP
delta_mfi_plot_log <- diff_inis %>% 
  filter(timepoint<=96 & plasmid=="SP505") %>% 
  filter(cellstate!="2i-LIFw") %>% 
  ggplot()+
  geom_point(aes(x=timepoint, y=(mfi_gfp_diff_log),color=factor(c1,levels=c("0ngDox","200ngDox",
  "400ngDox", "2000ngDox"))),pch=16,size = 1, alpha = 2/3)+
  geom_line(data=mean_diff_inis %>% filter(timepoint<=96&cellstate!="2i-LIFw"& plasmid=="SP505"),aes(x=timepoint, y=(mfi_gfp_diff_log_rep),
  color=factor(c1,levels=c("0ngDox","200ngDox","400ngDox", "2000ngDox"))))+  scale_fill_manual(values = brewer.pal(4, "Blues")) +
  scale_color_manual(values = brewer.pal(4, "Blues")) +
  scale_x_continuous(breaks=c(0,8,24,48,72,96), labels=c(0,0.33,1,2,3,4), expand = c(0, 0), limits=c(-4,100))+
  facet_wrap(~cellstate,ncol=4)+
  ylab("Delta log10(MFI GFP)")
fix1 <- set_panel_size(delta_mfi_plot_log, height = unit(2, "cm"), width = unit(1.5, "cm"))
ggsave(path = output_dir, "Fig5de_Delta_MFI_GFP.pdf", fix1, dpi = 300,
       useDingbats=FALSE)

# Delta MFI tdt
delta_mfi_plot_tdt_log <- diff_inis_tdt %>% 
  filter(timepoint<=96& plasmid=="SP505") %>% 
  filter(cellstate!="2i-LIFw") %>% 
  ggplot()+
  geom_point(aes(x=timepoint, y=(mfi_tdt_diff_log),color=factor(c1,levels=c("0ngDox","200ngDox","400ngDox", "2000ngDox"))),pch=16,size = 1, alpha = 2/3)+
  geom_line(data=mean_diff_inis_tdt %>% filter(timepoint<=96&cellstate!="2i-LIFw"& plasmid=="SP505"),aes(x=timepoint, y=(mfi_tdt_diff_log_rep),
  color=factor(c1,levels=c("0ngDox","200ngDox","400ngDox", "2000ngDox"))))+
  scale_fill_manual(values = brewer.pal(4, "Blues")) +
  scale_color_manual(values = brewer.pal(4, "Blues")) +
  scale_x_continuous(breaks=c(0,8,24,48,72,96), labels=c(0,0.33,1,2,3,4), expand = c(0, 0), limits=c(-4,100))+
  facet_wrap(~cellstate,ncol=4)+
  ylab("Delta log10(MFI tdT)")
fix1 <- set_panel_size(delta_mfi_plot_tdt_log, height = unit(2, "cm"), width = unit(1.5, "cm"))
ggsave(path = output_dir, "Fig5de_Delta_MFI_tdT.pdf", fix1, dpi = 300,
       useDingbats=FALSE)

