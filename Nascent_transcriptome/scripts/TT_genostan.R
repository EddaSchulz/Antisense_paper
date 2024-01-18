library(STAN)
library(tidyverse)
library(GenomicRanges)
library(egg)
library(gridExtra)
library(gplots)

args <- commandArgs(trailingOnly = TRUE)
day <- args[1]
data_dir <- args[2]

setwd(data_dir)

chr_vec <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

#Loads binned counts of TT-seq into R
binned_counts <- read.delim(paste0("XX_", day,"_binned_counts.txt"))  %>% 
  filter(chr %in% chr_vec) %>%
  mutate(chr = factor(chr, levels = chr_vec)) %>% 
  arrange(factor(chr, levels = chr_vec), start)
gRanges_binned <- GRanges(binned_counts[1:3])

#TTseq Matrix
TT_plus_matrix <- cbind(binned_counts[,4], binned_counts[,4])
colnames(TT_plus_matrix) <- c("TT_plus", "TT_plus")

TT_minus_matrix <- cbind(binned_counts[,5], binned_counts[,5])
colnames(TT_minus_matrix) <- c("TT_minus", "TT_minus")

TT_plus_list <- list(TT_plus_matrix)
names(TT_plus_list)  <- "TT_plus"

TT_minus_list <- list(TT_minus_matrix)
names(TT_minus_list)  <- "TT_minus"


#Running GenoStan model for plus and minus
nStates = 7
HMM_method = "PoissonLogNormal"

hmm.TT.plus = initHMM(TT_plus_list, nStates, HMM_method)
hmm.fitted.TT.plus = fitHMM(TT_plus_list, hmm.TT.plus, maxIters = 200, nCores = 1, verbose = TRUE)
viterbi_TT.plus = getViterbi(hmm.fitted.TT.plus, TT_plus_list)
avg_cov_TT.plus = getAvgSignal(viterbi_TT.plus, TT_plus_list)


hmm.TT.minus = initHMM(TT_minus_list, nStates, "PoissonLogNormal")
hmm.fitted.TT.minus = fitHMM(TT_minus_list, hmm.TT.minus, maxIters = 200, nCores = 1, verbose = TRUE)
viterbi_TT.minus = getViterbi(hmm.fitted.TT.minus, TT_minus_list)
avg_cov_TT.minus = getAvgSignal(viterbi_TT.minus, TT_minus_list)

## define state order and colors
heat = c("dark blue", "dodgerblue4", "darkred", "red", "orange", "gold", "yellow")
colfct = colorRampPalette(heat)
colpal_statemeans = colfct(200)

ord_TT.plus = order(apply(avg_cov_TT.plus,1,max), decreasing=TRUE)

statecols_TT.plus = rainbow(nStates)
names(statecols_TT.plus) = ord_TT.plus

pdf(paste0("Genostan_XX_", day, "_plus_heatmap.pdf"), 
    useDingbats = FALSE, onefile = FALSE)
heatmap.2(log(avg_cov_TT.plus+1)[as.character(ord_TT.plus),], margins=c(8,7), srtCol=45,
          RowSideColors=statecols_TT.plus[as.character(ord_TT.plus)], dendrogram="none",
          Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none",
          cellnote=round(avg_cov_TT.plus,1)[as.character(ord_TT.plus),], notecol="black")


dev.off()

ord_TT.minus = order(apply(avg_cov_TT.minus,1,max), decreasing=TRUE)

statecols_TT.minus = rainbow(nStates)
names(statecols_TT.minus) = ord_TT.minus

pdf(paste0("Genostan_XX_", day, "_minus_heatmap.pdf"), 
    useDingbats = FALSE, onefile = FALSE)
heatmap.2(log(avg_cov_TT.minus+1)[as.character(ord_TT.minus),], margins=c(8,7), srtCol=45,
          RowSideColors=statecols_TT.minus[as.character(ord_TT.minus)], dendrogram="none",
          Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none",
          cellnote=round(avg_cov_TT.minus,1)[as.character(ord_TT.minus),], notecol="black")


dev.off()

#Construct rgb bed file
factor_df_plus <- data.frame(statecols_TT.plus) %>% 
  rownames_to_column("state") 

factor_df_plus$rgb <- apply(factor_df_plus, 1, function(x) { 
  a <- as.vector(col2rgb(x[2]))
  str <- paste(a[1], a[2], a[3], sep = ",")
  return(str)
})

factor_df_plus_neat <- factor_df_plus %>% 
  transmute(names = state, itemRgb = rgb)

bed_df_plus <- data.frame(chr=seqnames(gRanges_binned),
                     start=start(gRanges_binned),
                     end=end(gRanges_binned),
                     names=viterbi_TT.plus$TT,
                     scores=c(rep("1000", length(gRanges_binned))),
                     strands=".",
                     thickStart = start(gRanges_binned), 
                     thickEnd = end(gRanges_binned)) %>% 
  left_join(factor_df_plus_neat) %>% 
  arrange(factor(chr, levels = chr_vec), start)

write.table(bed_df_plus, file=paste0("Genostan_XX_", day, "_plus_states.bed"), quote=F, 
            sep="\t", row.names=F, col.names=F)

factor_df_minus <- data.frame(statecols_TT.minus) %>% 
  rownames_to_column("state") 

factor_df_minus$rgb <- apply(factor_df_minus, 1, function(x) { 
  a <- as.vector(col2rgb(x[2]))
  str <- paste(a[1], a[2], a[3], sep = ",")
  return(str)
})

factor_df_minus_neat <- factor_df_minus %>% 
  transmute(names = state, itemRgb = rgb)

bed_df_minus <- data.frame(chr=seqnames(gRanges_binned),
                          start=start(gRanges_binned),
                          end=end(gRanges_binned),
                          names=viterbi_TT.minus$TT,
                          scores=c(rep("1000", length(gRanges_binned))),
                          strands=".",
                          thickStart = start(gRanges_binned), 
                          thickEnd = end(gRanges_binned)) %>% 
  left_join(factor_df_minus_neat) %>% 
  arrange(factor(chr, levels = chr_vec), start)

write.table(bed_df_minus, file=paste0("Genostan_XX_", day, "_minus_states.bed"), quote=F, 
            sep="\t", row.names=F, col.names=F)

#Reduce beds and construct rgb files
bed_df_reduced_plus <- bed_df_plus %>% 
  mutate(names = ifelse(names %in% ord_TT.plus[1:5], "plus", "none"), 
         itemRgb = ifelse(names == "plus", "255,128,0", "255,255,255"))

write.table(bed_df_reduced_plus, file=paste0("Genostan_XX_", day, "_plus_reduced.bed"), quote=F, 
            sep="\t", row.names=F, col.names=F)

bed_df_reduced_minus <- bed_df_minus %>% 
  mutate(names = ifelse(names %in% ord_TT.minus[1:5], "minus", "none"), 
         itemRgb = ifelse(names == "minus", "128,255,0", "255,255,255"))

write.table(bed_df_reduced_minus, file=paste0("Genostan_XX_", day, "_minus_reduced.bed"), quote=F, 
            sep="\t", row.names=F, col.names=F)

