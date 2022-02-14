library(ggplot2)
library(gridExtra)
library(MASS)
library(grDevices)
library(cowplot)
library(gridGraphics)


#############################
# 1. Confidence band plots
#############################

setwd("Enrichment-Inference-Supplemental-Materials\\simulation_files\\final_figure")
source("cp_bands_plots.r")
source("convenience.r")
load("cp_bands_all_dist.rdata")
cp_plot_res1 <- cp_band_plot(increment = c("Case 1", "Case 2", "Case 3", "Case 4", "Case 5"), axis.ylim = c(.15, .15), nox = 1)


#############################
# 2. Difference confidence band plots
#############################
source("cp_diff_bands_plots.r")

fn.v <- c("cp_diff_bands_binorm_nocorr.rdata", "cp_diff_bands_binorm_corr.rdata", "cp_diff_bands_bibeta_nocorr.rdata", "cp_diff_bands_bibeta_corr.rdata")
m.ls <- c()
k.true.ls <- list()
dist.v <- c("binorm", "binorm", "bibeta", "bibeta")
for(idx in seq_along(fn.v)) {
  print(fn.v[idx])
  
  load(fn.v[idx])
  
  m = 150000
  pi.0.true = 1/500
  params <- k.ls[[1]][[1]]$params
  r <- k.ls[[1]][[1]]$r
  dist <- dist.v[idx]
  k.diff.true <- computeTrueKDiff()
  
  k.true.ls <- c(k.true.ls, rep(list(k.diff.true), 6))
  
  tmp.ls <- c(k.ls[[1]][1:3], k.c.ls[[1]][1:3])
  tmp.ls <- tmp.ls[c(1,4,2,5,3,6)]
  m.ls <- c(m.ls, tmp.ls)
}

cp_plot_res2 <- cp_band_plot(m.ls, increment = c("Case 1", "Case 2", "Case 3", "Case 4"), axis.ylim = c(.15, .15), nox = 1)


pdf("simulation_plot2.pdf", height = 9, width = 7)
plot_grid(cp_plot_res1[[1]], cp_plot_res2[[1]], cp_plot_res1[[2]],  ncol = 1, rel_heights = c( .9, .9, .2))
dev.off()

