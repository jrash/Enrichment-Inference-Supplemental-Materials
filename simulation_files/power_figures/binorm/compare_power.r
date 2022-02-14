add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

allMetric <- function(ls, legend = 1) {
  simreps = ncol(ls[[1]][[1]]$covout)
  
  cols <- c("darkgreen","red", "purple", "darkorange", "cyan", "darkgrey", "black")
  frac <- log(ls[[1]][[1]]$r * 150000, base = 2)
  xlim <- c(min(frac), max(frac))
  mets <- 4
  
  ### Power
  
  plot(1, type="n", xaxt = 'n', ylab = "Power", xlab = "Compounds Tested",
       ylim = c(0,1), xlim = xlim)
  ticks <- c(1, 2, log(10, base = 2), seq(6, 8, by = 2), log(1500, base = 2), log(15000, base = 2))
  axis(1, at = ticks, labels = 2^ticks)
  
  for(i in 1:mets) {
    covprobs <- apply(ls[[1]][[i]]$covout, MARGIN=1, mean)
    
    # col.i <- ifelse(j == 2 & i == 4, cols[i+1], cols[i])
    lines(frac, covprobs, type="l", col = cols[i], lty = 1, lwd = 1.25)
    se <- sqrt(covprobs*(1-covprobs)/simreps)
    ucl <- covprobs + 1.96 * se
    lcl <- covprobs - 1.96 * se
    polygon(c(frac, rev(frac)), c(ucl, rev(lcl)),
            col = add.alpha(cols[i], .25), border = NA)
  }
  
  
  if(legend){
    legend(x= 0.5, y = 1.05, legend= c("EmpProc", "CorrBinom", "JZ Ind", "McNemar"),
           col=cols[1:4], lty=rep(1, 5), cex=.75, bty = "n", lwd = 2)
  }
  # par(old.par)
}