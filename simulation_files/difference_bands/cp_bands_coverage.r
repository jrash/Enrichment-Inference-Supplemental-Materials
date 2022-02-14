library(ggplot2)
library(gridExtra)
library(MASS)
library(grDevices)
library(cowplot)
setwd("G:\\My Drive\\JHO_lab\\recall_curve_papers\\JASA\\2021-09-07\\difference_cb")
source("cp_bands_plots.r")

######Convenience functions


##############################################

qmix <- function (p, params, pro = c(.5, .5), expand=.1, dist = "binorm") {
  nr <- 1000000
  XX <- rbinom(n=nr, size=1, prob=pro)
  if(dist == "binorm") {
    sd.s <- 1
    x = rnorm(n=nr, params[1], sd.s)*XX + rnorm(n=nr, params[2], sd.s)*(1-XX)
  } else if(dist == "bibeta") {
    x = rbeta(n=nr, params[[1]][1], params[[1]][2])*XX + rbeta(n=nr, params[[2]][1], params[[2]][2])*(1-XX)
  } 
  span <- seq(min(x) - expand * diff(range(x)), max(x) + expand * diff(range(x)), length = nr)
  cdf <- vector(mode = "numeric", length = nr)
  if(dist == "binorm") {
    sd.s <- 1
    cdf = pnorm(span, params[1], sd.s)*pro[1] + pnorm(span, params[2], sd.s)*(1-pro[1])
  } else if(dist == "bibeta") {
    cdf = pbeta(span, params[[1]][1], params[[1]][2])*pro[1] + pbeta(span, params[[2]][1], params[[2]][2])*(1-pro[1])
  } 
  quants <- stats::spline(cdf, span, method = "hyman", xout = p)$y
  quants[which(p < 0L | p > 1L)] <- NaN
  quants[which(p == 0L)] <- -Inf
  quants[which(p == 1L)] <- Inf
  if (any(is.nan(quants)))
    warning("Some quantile values could not be calculated. If all 'p's are within [0,1], try reducing the value of 'expand' and try again.")
  return(as.vector(quants))
}


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

computeTrueKDiff <- function(){
  
  if(dist == "binorm" | dist == "bibeta") {
    c1 <- qmix(1-r, params[[1]], c(pi.0.true, 1-pi.0.true), expand=.1, dist = dist)
    c2 <- qmix(1-r, params[[2]], c(pi.0.true, 1-pi.0.true), expand=.1, dist = dist)
  } else {
    stop("Wrong dist")
  }
  
  if (dist == "binorm") {
    sd.s = 1
    F.plus1 <- pnorm(c1, mean = params[[1]][1], sd.s)
    F.plus2 <- pnorm(c2, mean = params[[2]][1], sd.s)
  } else if (dist == "bibeta") {
    F.plus1 <- pbeta(c1, params[[1]][[1]][1], params[[1]][[1]][2])
    F.plus2 <- pbeta(c2, params[[2]][[1]][1], params[[2]][[1]][2])
  } 
  
  k.true1 <- vector(length = length(F.plus1))
  pi.true1 <- vector(length = length(F.plus1))
  k.true2 <- vector(length = length(F.plus2))
  pi.true2 <- vector(length = length(F.plus2))
  
  for (i in 1:length(k.true1)) {
    k.true1[i] <- 1 - F.plus1[i]
    pi.true1[i] <- k.true1[i]*pi.0.true/r[i]
    k.true2[i] <- 1 - F.plus2[i]
    pi.true2[i] <- k.true2[i]*pi.0.true/r[i]
  }
  
  # if k or pi is above 1, then they should actually be slightly below 1
  k.true1 <- ifelse(k.true1 > 1, 1 - 1E-10, k.true1)
  pi.true1 <- ifelse(pi.true1 > 1, 1 - 1E-10, pi.true1)
  k.true2 <- ifelse(k.true2 > 1, 1 - 1E-10, k.true2)
  pi.true2 <- ifelse(pi.true2 > 1, 1 - 1E-10, pi.true2)
  
  
  k.true1 - k.true2
}
##################################


fn.v <- c("power_band_binorm_nocorr.rdata", "power_band_binorm_corr.rdata", "power_band_bibeta_nocorr.rdata", "power_band_bibeta_corr.rdata")
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

pdf("cp_band_coverage.pdf", height = 9, width = 7)
cp_band_plot(m.ls, increment = c("Case 1", "Case 2", "Case 3", "Case 4"), axis.ylim = c(.15, .15), nox = 1)
dev.off()


