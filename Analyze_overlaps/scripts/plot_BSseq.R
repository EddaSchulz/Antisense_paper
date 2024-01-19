#Script to plot DNA methylation heatmaps
library(bsseq)
library(genomation)
library(egg)
library(gridExtra)
library(tidyverse)
library(viridis)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(),
                  axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
                  axis.text = element_text(size = 6), strip.text = element_text(size = 6)))


args <- commandArgs(trailingOnly = TRUE)
bs_dir <- args[1]
data_dir <- args[2]
fig_dir <- args[3]


setwd(bs_dir)
file_list = normalizePath(list.files(pattern="XX*.bedGraph$", full.names = TRUE))

XX_d0 = read.bismark("BSseq_XX_d0_CpG.bedGraph")
XX_d2 = read.bismark("BSseq_XX_d2_CpG.bedGraph")
XX_d4 = read.bismark("BSseq_XX_d4_CpG.bedGraph")

samples = bsseq::combine(XX_d0, XX_d2, XX_d4)

#Read Bed
setwd(data_dir)
XX_3prime_case <- GenomicRanges::promoters(readBed(paste0("XX_d0_3prime_case.bed")), 
                                              upstream = 500, downstream = 200)
XX_3prime_controls <- GenomicRanges::promoters(readBed(paste0("XX_d0_3prime_control.bed")), 
                                                  upstream = 500, downstream = 200)
XX_5prime_case <- GenomicRanges::promoters(readBed(paste0("XX_d0_5prime_case.bed")), 
                                              upstream = 500, downstream = 200)
XX_5prime_controls <- GenomicRanges::promoters(readBed(paste0("XX_d0_5prime_control.bed")),
                                                  upstream = 500, downstream = 200)


#Associate regions with meth%
meth_3prime_case <- getMeth(samples, XX_3prime_case, type = "raw", what = "perRegion")
colnames(meth_3prime_case) <- c("d0", "d2", "d4")
meth_3prime_controls <- getMeth(samples, XX_3prime_controls, type = "raw", what = "perRegion")
colnames(meth_3prime_controls) <- c("d0", "d2", "d4")

meth_5prime_case <- getMeth(samples, XX_5prime_case, type = "raw", what = "perRegion")
colnames(meth_5prime_case) <- c("d0", "d2", "d4")
meth_5prime_controls <- getMeth(samples, XX_5prime_controls, type = "raw", what = "perRegion")
colnames(meth_5prime_controls) <- c("d0", "d2", "d4")

#Turn matrices into dataframes
meth_3prime_case <- data.frame(meth_3prime_case) %>% 
  mutate(location = "3prime")
meth_3prime_controls<- data.frame(meth_3prime_controls) %>% 
  mutate(location = "controls")

meth_3prime <- rbind(meth_3prime_case, meth_3prime_controls) %>% 
  mutate(group = "3prime")


meth_5prime_case <- data.frame(meth_5prime_case) %>% 
  mutate(location = "5prime")
meth_5prime_controls<- data.frame(meth_5prime_controls) %>% 
  mutate(location = "controls")

meth_5prime <- rbind(meth_5prime_case, meth_5prime_controls) %>% 
  mutate(group = "5prime")

meth_total <- rbind(meth_3prime, meth_5prime)

#Plots methylation% over time 
setwd(fig_dir)

meth_viol_3prime <- meth_total %>%
  filter(group == "3prime") %>% 
  pivot_longer(c(1:3), names_to = "day", values_to = "meth") %>% 
  ggplot(aes(y = meth, x = day, fill = location)) +
  geom_split_violin(scale = "width") +
  scale_fill_manual(values = c("#402F87", "#777777"))

fix <- set_panel_size(meth_viol_3prime, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("Fig4F_d0_matched_meth_violin_3prime.pdf", fix, dpi = 300,
       useDingbats=FALSE)

meth_viol_5prime <- meth_total %>%
  filter(group == "5prime") %>% 
  pivot_longer(c(1:3), names_to = "day", values_to = "meth") %>% 
  ggplot(aes(y = meth, x = day, fill = location)) +
  geom_split_violin(scale = "width") +
  scale_fill_manual(values = c("#F7A539", "#777777"))

fix <- set_panel_size(meth_viol_5prime, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("Fig4F_d0_matched_meth_violin_5prime.pdf", fix, dpi = 300,
       useDingbats=FALSE)
