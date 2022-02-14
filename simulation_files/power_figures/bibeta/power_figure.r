setwd("G:\\My Drive\\JHO_lab\\recall_curve_papers\\JASA\\2021-09-07\\power_figures\\bibeta")
source("compare_power.r")
source("convenience.r")

pdf("power_figure.pdf")
pi.0 = 1/500; m = 150000; rlim = .1;
plus2 = F; lift.ylim = 50;
distr = "bibeta";
par(mfrow = c(2, 2), mar = c(4, 4, 3, 3))
set.seed(125)

col.list <- list(c("blue", "orange"), c("green", "orange"))
params.ls <- list(list(c(5, 2), c(2, 5)), list(c(4, 2), c(2, 5)))

c <- seq(0, 1, .001)

F.plus1 <- pbeta(c, params.ls[[1]][[1]][1], params.ls[[1]][[1]][2])
F.plus2 <- pbeta(c, params.ls[[2]][[1]][1], params.ls[[2]][[1]][2])
F.neg <- pbeta(c, params.ls[[1]][[2]][1], params.ls[[1]][[2]][2])
f.plus1 <- dbeta(c, params.ls[[1]][[1]][1], params.ls[[1]][[1]][2])
f.plus2 <- dbeta(c, params.ls[[2]][[1]][1], params.ls[[2]][[1]][2])
f.neg <- dbeta(c, params.ls[[1]][[2]][1], params.ls[[1]][[2]][2])

# PDFs

plot(c, f.plus1, xlab="s", type = "l",
     ylab="f(s)", col = "blue", cex.lab=1.3)
lines(c, f.plus2, xlab="s", type = "l",
      ylab="f(s)", col = "orange")
lines(c, f.neg, xlab="s", type = "l",
      ylab="f(s)", col = "black")
legend(x="top", legend=c("Algorithm 1","Algorithm 2"), 
       col=c(col.list[[1]][1], col.list[[1]][2]),
       lty=c(1,1), cex=.75, bty = "n", lwd = 2)
fig_label("A", cex=2.5, region="figure", pos="topleft")

xlim = log2(c(1/150000, .1)*150000)
for(metric in  c("k")) {
  for(pi.0 in c(1/500)){
    if (metric == "k") {
      plot(1, type="n", xaxt = 'n', ylab = "Recall", xlab = "Expected Compounds Tested (n = 150,000)",
           ylim = c(0,1), xlim = xlim)
      ticks <- c(1, 2, log(10, base = 2), seq(6, 8, by = 2), log(1500, base = 2), log(15000, base = 2))
      axis(1, at = ticks, labels = 2^ticks)
      fig_label("B", cex=2.5, region="figure", pos="topleft")
    } 
    for(i in 1:2) {
      params <- params.ls[[i]]
      #print(distr)
      if (distr == "binorm") {
        c <- seq(-10, 10, .0001)
        F.plus <- pnorm(c, mean = params[1], sd = 1)
        F.neg <- pnorm(c, mean = params[2], sd = 1)
        f.plus <- dnorm(c, mean = params[1], sd = 1)
        f.neg <- dnorm(c, mean = params[2], sd = 1)
      } else if (distr == "bibeta") {
        c <- seq(0, 1, .001)
        F.plus <- pbeta(c, params[[1]][1], params[[1]][2])
        F.neg <- pbeta(c, params[[2]][1], params[[2]][2])
        f.plus <- dbeta(c, params[[1]][1], params[[1]][2])
        f.neg <- dbeta(c, params[[2]][1], params[[2]][2])
      } else if (distr == "unif") {
        c <- seq(0, 1, .00001)
        F.plus <- punif(c, params[[1]][1], params[[1]][2])
        F.neg <- punif(c, params[[2]][1], params[[2]][2])
        f.plus <- dunif(c, params[[1]][1], params[[1]][2])
        f.neg <- dunif(c, params[[2]][1], params[[2]][2])
      }
      
      r <- 1 - pi.0 * F.plus - (1 - pi.0) * F.neg
      # only look at beginning of the curve
      idx <- r < rlim & (r > 1/m)
      r <- r[idx]
      frac <- log(150000*r, base = 2)
      k <- 1 - F.plus[idx]
      pi <- k*pi.0/r
      if (plus2){
        # Expected number of actives tested divided by expected number of actives
        k_c <- (m*pi.0*k+2)/(m*pi.0+4)
        # Expected number of actives tested divided by expected number of tested
        pi_c <- (m*r*pi+2)/(m*r+4)
      }
      f.s <- pi.0*f.plus[idx] + (1-pi.0)*f.neg[idx]
      lam <- pi.0*f.plus[idx]/f.s
      
      
      if (plus2){
        var.k.ratio <- 1-2*lam+(lam^2*(1-r))/(pi_c*(1-k_c))
        var.k.b <- ((m*pi.0)^-1)*k_c*(1-k_c)
        var.k <- var.k.b*var.k.ratio
      } else {
        var.k.ratio <- 1-2*lam+(lam^2*(1-r))/(pi*(1-k))
        var.k.b <- ((m*pi.0)^-1)*k*(1-k)
        var.k <- var.k.b*var.k.ratio
      }
      
      if (metric == "k") {
        lines(frac, k, lty = 1, lwd=1, col = col.list[[1]][i])
        ucl <- k + 1.96*sqrt(var.k); ucl <- ifelse(ucl > 1, 1, ucl)
        lcl <- k - 1.96*sqrt(var.k); lcl <- ifelse(lcl < 0, 0, lcl)
        polygon(c(frac, rev(frac)), c(ucl, rev(lcl)),
                col = add.alpha(col.list[[1]][i], .25), border = NA)
        lines(frac, 2^frac/150000, col = "grey")
      }
    }
  }
}

fixK <- function(ls1, ls2){
  for(i in 1:3) {
    ls1[[i]][[1]] <- ls2[[i]][[1]]
  }
  ls1
}


load("power_all_bibeta_unpooled_nocorr.rdata")
ls <- k.ls

allMetric(ls, legend = 1)
fig_label("C", cex=2.5, region="figure", pos="topleft")

load("power_all_bibeta_unpooled_corr.rdata")
ls <- k.ls

allMetric(ls, legend = 0)
fig_label("D", cex=2.5, region="figure", pos="topleft")

dev.off()
