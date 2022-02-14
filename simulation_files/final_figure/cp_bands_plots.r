cp_band_plot <- function(increment, metric, axis.ylim = c(.15, .15), nox = 0) {
  metric = "k"
  if (metric == "k") {
    m.ls <- k.ls
  } else if (metric == "pi") {
    m.ls <- pi.ls
  }
  covprobs <- vector(length = length(m.ls))
  for(i in seq_along(m.ls)) {
    covprobs[i] <- mean(m.ls[[i]]$covout)
  }
  simreps = length(m.ls[[1]]$interval.ls)
  se <- sqrt(covprobs*(1-covprobs)/simreps)
  ucl <- covprobs + 1.96 * se
  lcl <- covprobs - 1.96 * se
  increment <- factor(increment, levels = increment)
  method.v <- c("Sup-t", "Sup-t Plus", "Bonferroni", "Bonferroni Plus",
                "Point-Wise", "Point-Wise Plus")
  plot.df <- data.frame(covprobs, ucl, lcl, method.v, rep(increment, each = length(method.v)))
  colnames(plot.df) <- c("CP", "UCL", "LCL", "Method", "Increment")
  # plot.df <- plot.df[!plot.df$Method %in% c("Sup-t", "Bonferroni"), ]
  
  
  pd <- position_dodge(0.2)
  k<-ggplot(plot.df,aes(x=Increment,y=CP, group=Method,color=Method))
  k1 <- k+ geom_errorbar(aes(ymin=LCL, ymax=UCL), width=1, position=pd, size = .5) +
    geom_point(position=pd, size = 3, shape = 23) +
    ylab("Coverage Probability") +
    xlab(label = switch(nox, NULL, "Number of Grid Points")) +
    scale_y_continuous(breaks=seq(0, 1, .2)) +         # Set tick every 4
    theme_bw() +
    theme(legend.position = "none") + 
    geom_hline(yintercept=.95, linetype = "dashed") + 
    scale_colour_brewer(type = "qual", palette = "Paired")
  
  
  expwidth <- vector(length = length(m.ls))
  ucl.p <- vector(length = length(m.ls))
  lcl.p <- vector(length = length(m.ls))
  for(i in seq_along(m.ls)) {
    for(j in 1:simreps) {
      lcl <- m.ls[[i]]$interval.ls[[j]][, 1]
      lcl <- ifelse(lcl < 0, 0, lcl)
      ucl <- m.ls[[i]]$interval.ls[[j]][, 2]
      ucl <- ifelse(ucl > 1, 1, ucl)
      m.ls[[i]]$interval.ls[[j]][, 1] <- lcl
      m.ls[[i]]$interval.ls[[j]][, 2] <- ucl
      m.ls[[i]]$widthout[, j] <- ucl - lcl
    }
    expwidth[i] <- mean(m.ls[[i]]$widthout)
    se <- sd(m.ls[[i]]$widthout)/sqrt(simreps*length(increment))
    ucl.p[i] <- expwidth[i] + 1.96 * se
    lcl.p[i] <- expwidth[i] - 1.96 * se
  }
  increment <- factor(increment, levels = increment)
  plot.df <- data.frame(expwidth, ucl.p, lcl.p, as.character(method.v), rep(increment, each = length(method.v)))
  colnames(plot.df) <- c("EW", "UCL", "LCL", "Method", "Increment")
  # plot.df <- plot.df[!plot.df$Method %in% c("Sup-t", "Bonferroni"), ]
  
  pd <- position_dodge(0.2)
  k<-ggplot(plot.df,aes(x=Increment,y=EW, group=Method,color=Method))
  
  k2 <- k+ geom_errorbar(aes(ymin=LCL, ymax=UCL), width=1, position=pd, size = .5) +
    geom_point(position=pd, size = 3, shape = 23) +
    ylab("Average Width") +
    xlab(label = switch(nox, NULL, "Number of Grid Points")) +
    coord_cartesian(ylim = c(0, axis.ylim[1])) +
    scale_y_continuous(breaks = seq(0, axis.ylim[1], length.out = 4)) + 
    theme_bw() +
    scale_fill_discrete(name = "Methods", labels = method.v) +
    theme(legend.position="none") + 
    scale_colour_brewer(type = "qual", palette = "Paired")
  
  prow <- plot_grid(
    k1 + theme(legend.position="none"),
    k2 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    nrow = 1,
    ncol = 2
  )
  legend_b <- get_legend(
    k1 + 
      guides(color = guide_legend(nrow = 2)) +
      theme(legend.position = "bottom", legend.title=element_text(size=10), 
            legend.text=element_text(size=8))
  )
  
  return(list(prow, legend_b))
  
}