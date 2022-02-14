library(ggplot2)
library(gridExtra)
library(MASS)
library(grDevices)
library(cowplot)
library(gridGraphics)
library(RColorBrewer)

setwd("Enrichment-Inference-Supplemental-Materials\\simulation_files\\final_figure")
source("convenience.r")

############
# 1. Power plots
############

load("power_all_binorm_unpooled_corr.rdata")
k1.ls <- k.ls
k2.ls <- k.c.ls

m = 150000
pi.0.true = 1/500
params <- k1.ls[[1]][[1]]$params
r <- k1.ls[[1]][[1]]$r
dist <- "binorm"

k.diff.true <- computeTrueKDiff()

source("compare_power.r")
source("convenience.r")

powerplot <- function(){
  par(mfrow = c(1, 2), mai = c(1.02, 0.82, 0.25, 0.2), oma = c(0, 0, 1, 0))
  load("power_all_bibeta_unpooled_nocorr.rdata")
  allMetric(ls = k.ls, legend = 1)
  fig_label("A", cex=1.3, region="figure", pos="topleft", font = 2)
  load("power_all_bibeta_unpooled_corr.rdata")
  allMetric(ls = k.ls, legend = 0)
  fig_label("B", cex=1.3, region="figure", pos="topleft", font = 2)
}


#############################
# 2. Coverage Probability Plots
#############################

covProbPlot <- function() {
  
  par(mfrow = c(1, 2), mai = c(1.02, 0.82, 0.25, 0.2), oma = c(0, 0, 1, 0))
  
  all.k <- list(k1.ls, k2.ls)
  
  cols <- rep(brewer.pal(8,"Paired")[c(2, 4, 6, 8)], each = 2)
  type.v <- rep(c(2, 1), 4)
  frac <- log(all.k[[1]][[1]][[1]]$r * 150000, base = 2)
  xlim <- c(min(frac), max(frac))
  
  plot(1, type="n", xaxt = 'n', ylab = "Coverage Probability", xlab = "Compounds Tested",
       ylim = c(.5,1), xlim = xlim, main = "", cex.lab = 1)
  ticks <- c(1, 2, log(10, base = 2), seq(6, 8, by = 2), log(1500, base = 2), log(15000, base = 2))
  axis(1, at = ticks, labels = 2^ticks)
  
  simrep <- length(all.k[[1]][[1]][[1]]$interval.ls)
  idx <- 1
  for(met in c(1, 4, 3)) {
    for(i in 1:2){
      cov.out <- vector(mode = "numeric", length = length(k.diff.true))
      for(l in 1:simrep){
        CI.int <- (all.k[[i]][[1]][[met]]$interval.ls)[[l]]
        for(j in seq_along(k.diff.true)) {
          cov.out[j] <- cov.out[j] + ((round(CI.int[j, 1], 10) <= k.diff.true[j]) & (round(CI.int[j, 2], 10) >= k.diff.true[j]))
        }
      }
      covprobs <- cov.out/simrep

      lines(frac, covprobs, type="l", col = cols[idx], lty = type.v[idx], lwd = 2)
      se <- sqrt(covprobs*(1-covprobs)/simrep)
      ucl <- covprobs + 1.96 * se
      lcl <- covprobs - 1.96 * se
      polygon(c(frac, rev(frac)), c(ucl, rev(lcl)),
              col = add.alpha(cols[idx], .25), border = NA)
      idx <- idx + 1
    }
  }
  abline(h = .95, lty = "dashed")
  legend(x = frac[length(frac) - 13], y = .725, legend=c("EmProc", "CorrBinom","IndJZ"),
         col=cols[c(1, 3, 5)], lty=rep(1, 3), cex = .75, lwd = rep(2, 3), bty = "n")
  legend(x = frac[length(frac) - 13], y = .625, legend = c("Unadjusted", "Plus Adjusted"), col = c("black", "black"), lty = c(2, 1), cex = .75, bty = "n", lwd = 2)
  fig_label("C", cex=1.3, region="figure", pos="topleft", font = 2)
  
  plot(1, type="n", xaxt = 'n', ylab = "Average Width", xlab = "Compounds Tested",
       ylim = c(0, .1), xlim = xlim, cex.lab = 1)
  ticks <- c(1, 2, log(10, base = 2), seq(6, 8, by = 2), log(1500, base = 2), log(15000, base = 2))
  axis(1, at = ticks, labels = 2^ticks, cex.axis = 1)
  
  idx <- 1
  for (met in c(1, 4, 3)) {
    for(i in 1:2){
      widthout <- matrix(ncol = length(k.diff.true), nrow = simrep)
      for(l in 1:simrep){
        CI.int <- (all.k[[i]][[1]][[met]]$interval.ls)[[l]]
        widthout[l, ] <- (CI.int[, 2] - CI.int[, 1])
      }
      expwidth <- apply(widthout, MARGIN=2, mean)
      se <- apply(widthout, MARGIN=2, sd)/sqrt(simrep)
      lines(frac, expwidth, type="l", col = cols[idx], lty = type.v[idx], lwd = 1.5)
      ucl <- expwidth + 1.96 * se
      lcl <- expwidth - 1.96 * se
      polygon(c(frac, rev(frac)), c(ucl, rev(lcl)),
              col = add.alpha(cols[idx], .25), border = NA)
      idx <- idx + 1
    }
  }

  fig_label("D", cex=1.3, region="figure", pos="topleft", font = 2)
  
}


#############################
# 3. Make Plot
#############################

pdf("simulation_plot1.pdf", height = 9, width = 7)

plot_grid(powerplot, covProbPlot, ncol = 1, rel_heights = c( 1, 1))

dev.off()
