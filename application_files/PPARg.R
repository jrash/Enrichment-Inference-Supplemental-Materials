## This file creates all figures, and the table, for the PPARg application
## Figure 1: hit enrichment curves, without confidence bands
## Figure 5: hit enrichment curves, with confidence bands
## Figure 6: hit enrichment differences, with confidence bands
## Table 2: pointwise hypothesis testing results
## Figure 4: pointwise hypothesis testing & confidence intervals

# Get most up-to-date version of chemmodlab from github ....
# detach("package:chemmodlab", force = TRUE)
# devtools::install_github("jrash/chemmodlab")
library("chemmodlab")
data(pparg)

## Figure 1: hit enrichment curves for PPARg application ----
pdf("Fig1.pdf", height = 4, width = 7)
HitEnrich(S.df = pparg[,c(2,5,8,11,14)], y = pparg[,3], x.max = NULL, labels = 
            c("Surflex-dock", "ICM","Vina", "Minimum Rank", "Maximum z-score"), 
          log = T)
dev.off()

## Figure 5: confidence bands for PPARg ----
pdf("Fig5.pdf", height = 4, width = 7)
HitEnrich(S.df = pparg[,c(14,2,5)], y = pparg[,3], x.max = NULL, labels = 
            c("Maximum z-score", "Surflex-dock", "ICM"), 
          log = T, conf = T, conf.level = .95)
dev.off()

## Figure 6: confidence bands of differences for PPARg ----
p.ls <- list()
p.ls[[1]] <- ~{
  HitEnrichDiff(S.df = pparg[,c(14,2)], y = pparg[,3], x.max = NULL, labels = 
                  c("Maximum z-score","Surflex-dock"), 
                log = T, conf.level = .95)
}
p.ls[[2]] <- ~{
  HitEnrichDiff(S.df = pparg[,c(14,5)], y = pparg[,3], x.max = NULL, labels = 
                  c("Maximum z-score","ICM"), 
                log = T, conf.level = .95)
}
p.ls[[3]] <- ~{
  HitEnrichDiff(S.df = pparg[,c(2,5)], y = pparg[,3], x.max = NULL, labels = 
                  c("Surflex-dock","ICM"), 
                log = T, conf.level = .95)
}
library("cowplot")
prow <- plot_grid(p.ls[[1]], NULL, p.ls[[2]], p.ls[[3]], ncol = 2,
                  labels = c("A", "", "B", "C"))

pdf("Fig6.pdf", height = 4, width = 6, pointsize = 9)
prow
dev.off()



## Table 2: pointwise hypothesis testing results ----
n <- nrow(pparg)
tested <- c(3,32,321)
ntested <- length(tested)
dock <- c("maxz","surf","icm"); ndock <- length(dock)
mthd <- c("EmProc","binomial","mcnemar","JZ ind"); nmthd <- length(mthd)
Rdiff <- array( dim = c(ndock, ntested, ndock, nmthd),
                dimnames = list( dock, as.character(tested), 
                                 dock, mthd )  )
RdiffSE <- array( dim = c(ndock, ntested, ndock, nmthd),
                  dimnames = list( dock, as.character(tested), 
                                   dock, mthd )  )
pvals <- array( dim = c(ndock, ntested, ndock, nmthd),
                dimnames = list( dock, as.character(tested), 
                                 dock, mthd )  )
cilow <- array( dim = c(ndock, ntested, ndock, nmthd),
                dimnames = list( dock, as.character(tested), 
                                 dock, mthd )  )
cihigh <- array( dim = c(ndock, ntested, ndock, nmthd),
                 dimnames = list( dock, as.character(tested), 
                                  dock, mthd )  )



# pvalue etc
for(i in 1:(ndock-1)){
  for(j in (i+1):ndock){
    for(k in 1:nmthd){
      keepit <- 
        PerfCurveTest( S1=eval(parse(text=paste("pparg$",dock[i],"_scores",sep=""))),
                       S2=eval(parse(text=paste("pparg$",dock[j],"_scores",sep=""))),
                       X=pparg$surf_actives, r=tested/n, metric="rec",
                       method=mthd[k] , plus=F, pool=F, alpha=.05)
      keepit2 <- 
        PerfCurveTest( S1=eval(parse(text=paste("pparg$",dock[i],"_scores",sep=""))),
                       S2=eval(parse(text=paste("pparg$",dock[j],"_scores",sep=""))),
                       X=pparg$surf_actives, r=tested/n, metric="rec",
                       method=mthd[k] , plus=T, pool=F, alpha=.05)
      Rdiff[i,,j,k] <- keepit$diff_estimate
      RdiffSE[i,,j,k] <- keepit$std_err
      pvals[i,,j,k] <- keepit$p_value
      cilow[i,,j,k] <- keepit2$ci_interval[,1]
      cihigh[i,,j,k] <- keepit2$ci_interval[,2]
    }
  }
}

# Restructure results into data frame ----
Rdiff_df <- data.frame(Rdiff[1,,2,])
Rdiff_df <- rbind(Rdiff_df, data.frame(Rdiff[1,,3,]))
Rdiff_df <- rbind(Rdiff_df, data.frame(Rdiff[2,,3,]))

RdiffSE_df <- data.frame(RdiffSE[1,,2,])
RdiffSE_df <- rbind(RdiffSE_df, data.frame(RdiffSE[1,,3,]))
RdiffSE_df <- rbind(RdiffSE_df, data.frame(RdiffSE[2,,3,]))

pval_df <- data.frame(pvals[1,,2,])
pval_df <- rbind(pval_df, data.frame(pvals[1,,3,]))
pval_df <- rbind(pval_df, data.frame(pvals[2,,3,]))

cilow_df <- data.frame(cilow[1,,2,])
cilow_df <- rbind(cilow_df, data.frame(cilow[1,,3,]))
cilow_df <- rbind(cilow_df, data.frame(cilow[2,,3,]))

cihigh_df <- data.frame(cihigh[1,,2,])
cihigh_df <- rbind(cihigh_df, data.frame(cihigh[1,,3,]))
cihigh_df <- rbind(cihigh_df, data.frame(cihigh[2,,3,]))


# Function movearound restructures the p-values and such to display in a table ----
movearound <- function( dfpval, dfCIl, dfCIu ){
  table_rec <- data.frame( 
    Comparison = c(rep("maxz - surf", 3), rep("maxz - icm", 3), rep("surf - icm", 3)),
    Difference = (dfCIu + dfCIl) / 2,
    PlusMinus = (dfCIu - dfCIl) / 2,
    Raw.p = dfpval, 
    Adj.p = p.adjust(as.vector(dfpval), method = "BH"),
    Fraction = rep(tested/n, 3) )
  
  table_EF <- table_rec
  table_EF$Difference <- table_EF$Difference/table_EF$Fraction
  table_EF$PlusMinus <- table_EF$PlusMinus/table_EF$Fraction
  
  table_rec <- table_rec %>% 
    pivot_wider( names_from = "Fraction", 
                 values_from = c("Difference", "PlusMinus", "Raw.p", "Adj.p") )
  ## The ordering of pivot_wider is not good and there is an open issue
  ## https://github.com/tidyverse/tidyr/issues/839
  table_rec <- table_rec[, c(1, c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12) + 1)]
  table_rec[, -1] <- signif(table_rec[, -1], 3)
  
  table_EF <- table_EF %>% 
    pivot_wider( names_from = "Fraction", 
                 values_from = c("Difference", "PlusMinus", "Raw.p", "Adj.p"))
  table_EF <- table_EF[, c(1, c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12) + 1)]
  table_EF[, -1] <- signif(table_EF[, -1], 3)
  
  #  return(list(table_rec, table_EF))
  return(table_rec)
}

table_rec_EmProc <- movearound( pval_df$EmProc, cilow_df$EmProc, cihigh_df$EmProc )
table_rec_mcnemar <- movearound( pval_df$mcnemar, cilow_df$mcnemar, cihigh_df$mcnemar )
table_rec_JZ.ind <- movearound( pval_df$JZ.ind, cilow_df$JZ.ind, cihigh_df$JZ.ind )
table_rec_binomial <- movearound( pval_df$binomial, cilow_df$binomial,
                                  cihigh_df$binomial )

movearound <- function( dfdiff, dfse, dfpval ){
  table_rec <- data.frame( 
    Comparison = c(rep("maxz - surf", 3), rep("maxz - icm", 3), rep("surf - icm", 3)),
    Difference = dfdiff,
    StdErr = dfse,
    Raw.p = dfpval, 
    Adj.p = p.adjust(as.vector(dfpval), method = "BH"),
    Fraction = rep(tested/n, 3) )
  
  table_EF <- table_rec
  table_EF$Difference <- table_EF$Difference/table_EF$Fraction
  table_EF$StdErr <- table_EF$StdErr/table_EF$Fraction
  
  table_rec <- table_rec %>% 
    pivot_wider( names_from = "Fraction", 
                 values_from = c("Difference", "StdErr", "Raw.p", "Adj.p") )
  ## The ordering of pivot_wider is not good and there is an open issue
  ## https://github.com/tidyverse/tidyr/issues/839
  table_rec <- table_rec[, c(1, c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12) + 1)]
  table_rec[, -1] <- signif(table_rec[, -1], 3)
  
  table_EF <- table_EF %>% 
    pivot_wider( names_from = "Fraction", 
                 values_from = c("Difference", "StdErr", "Raw.p", "Adj.p"))
  table_EF <- table_EF[, c(1, c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12) + 1)]
  table_EF[, -1] <- signif(table_EF[, -1], 3)
  
  #  return(list(table_rec, table_EF))
  return(table_rec)
}


table_rec_EmProc <- movearound( Rdiff_df$EmProc, RdiffSE_df$EmProc, pval_df$EmProc )
table_rec_mcnemar <- movearound( Rdiff_df$mcnemar, 
                                 RdiffSE_df$mcnemar, pval_df$mcnemar )
table_rec_JZ.ind <- movearound( Rdiff_df$JZ.ind, RdiffSE_df$JZ.ind, pval_df$JZ.ind )
table_rec_binomial <- movearound( Rdiff_df$binomial, RdiffSE_df$binomial, 
                                  pval_df$binomial )

# Create data frame of all methods ----
table_df <- rbind(apply(table_rec_EmProc, 2, as.character), 
                  apply(table_rec_mcnemar, 2, as.character),
                  apply(table_rec_JZ.ind, 2, as.character),
                  apply(table_rec_binomial, 2, as.character) )
colnames(table_df) <- c("", rep(c("Difference", "Std Err", "Raw p", "Adj p"), 3))


# Now actually create the table ----
library("tidyr")
library("kableExtra")
kable(table_df, format = "latex", align = c('l', rep('c',times=12)), booktabs = T,
      caption = "Pairwise comparison of scoring methods surflex-dock (surf), ICM, and their consensus (maxz), at three testing fractions decided a priori. For each pair of scoring methods, differences and standard errors of estimated hit enrichment curves are shown, along with raw and Benjamini-Hochberg adjusted $p$-values.\\label{tab:applic}", escape = T ) %>% 
  kable_styling(latex_options = c("scale_down"), position = "center") %>% 
  add_header_above(c(" ", "N tested = 3" = 4, "N tested = 32" = 4, 
                     "N tested = 321" = 4), italic = T, bold = T) %>% 
  pack_rows(index=c("EmProc" = 3, "McNemar" = 3, "IndJZ" = 3, "CorrBinom" = 3)) %>% 
  column_spec(c(5, 9), border_right = T) 



## Figure 4: pointwise hypothesis testing & confidence intervals ----
col.cd <- c("blue","red","cyan", "orange", "darkgreen")

# Calculate p-values and confidence intervals ----------------------
n <- nrow(pparg)
xrange <- log( c(.001,.1), base=2 )
xvals <- seq(xrange[1],xrange[2],length.out=15)
#tested <- ceiling( ( 2^xvals )* n )
tested <- unique( ceiling( ( 2^xvals )* n ) )
ntested <- length(tested)
dock <- c("maxz","surf","icm"); ndock <- length(dock)
mthd <- c("EmProc","binomial","mcnemar","JZ ind"); nmthd <- length(mthd)
recdiff <- array( dim = c(ndock, ntested, ndock, nmthd),
                  dimnames = list( dock, as.character(tested), 
                                   dock, mthd )  )
pvals <- array( dim = c(ndock, ntested, ndock, nmthd),
                dimnames = list( dock, as.character(tested), 
                                 dock, mthd )  )
cilow <- array( dim = c(ndock, ntested, ndock, nmthd),
                dimnames = list( dock, as.character(tested), 
                                 dock, mthd )  )
cihigh <- array( dim = c(ndock, ntested, ndock, nmthd),
                 dimnames = list( dock, as.character(tested), 
                                  dock, mthd )  )
for(i in 1:(ndock-1)){
  for(j in (i+1):ndock){
    for(k in 1:nmthd){
      keepit <- PerfCurveTest(
        S1=eval(parse(text=paste("pparg$",dock[i],"_scores",sep=""))),
        S2=eval(parse(text=paste("pparg$",dock[j],"_scores",sep=""))),
        X=pparg$surf_actives, r=tested/n, metric="rec",
        method=mthd[k] , plus=F, pool=F, alpha=.05)
      keepit2 <- 
        PerfCurveTest( S1=eval(parse(text=paste("pparg$",dock[i],"_scores",sep=""))),
                       S2=eval(parse(text=paste("pparg$",dock[j],"_scores",sep=""))),
                       X=pparg$surf_actives, r=tested/n, metric="rec",
                       method=mthd[k] , plus=T, pool=F, alpha=.05)
      recdiff[i,,j,k] <- keepit$diff_estimate
      pvals[i,,j,k] <- keepit$p_value
      cilow[i,,j,k] <- keepit2$ci_interval[,1]
      cihigh[i,,j,k] <- keepit2$ci_interval[,2]
    }
  }
}

# Intialize plot list
p.ls <- list()

# Density plots -----------------
library(ggpubr)
dens1 <- density(pparg$maxz_scores[pparg$surf_actives==1], from=-2.5,to=3.5)
dens0 <- density(pparg$maxz_scores[pparg$surf_actives==0], from=-2.5,to=3.5)
# install.packages("LaplacesDemon")
KLscore <- LaplacesDemon::KLD(dens0$y, dens1$y)$sum.KLD.py.px
p.ls[[1]] <- ggplot(pparg, aes(x = maxz_scores, fill = (surf_actives==1))) +
  geom_density(alpha = 0.5) + theme(legend.position="none") +
  annotate(geom="label", x=-2.3, y=0.6, label="Consensus", hjust = 0) +
  annotate(geom="text", x=-2.3, y=0.4, label=paste("KL=",signif(KLscore,dig=3)), 
           hjust = 0, size=2) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6),
        axis.text=element_text(size=6))

dens1 <- density(pparg$surf_scores[pparg$surf_actives==1], from=-40,to=20)
dens0 <- density(pparg$surf_scores[pparg$surf_actives==0], from=-40,to=20)
# KLscore <- LaplacesDemon::KLD(dens0$y, dens1$y)
# plot(KLscore$KLD.py.px)
KLscore <- LaplacesDemon::KLD(dens0$y, dens1$y)$sum.KLD.py.px
# KLscore
p.ls[[2]] <- ggplot(pparg, aes(x = surf_scores, fill = (surf_actives==1))) +
  geom_density(alpha = 0.5) + theme(legend.position="none") +
  annotate(geom="label", x=-40, y=0.175, label="Surflex-Dock", hjust = 0) +
  annotate(geom="text", x=-40, y=0.125, label=paste("KL=",signif(KLscore,dig=3)), 
           hjust = 0, size=2) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6),
        axis.text=element_text(size=6))

dens1 <- density(pparg$icm_scores[pparg$surf_actives==1], from=-350,to=100)
dens0 <- density(pparg$icm_scores[pparg$surf_actives==0], from=-350,to=100)
# KLscore <- LaplacesDemon::KLD(dens0$y, dens1$y)
# plot(KLscore$KLD.py.px)
KLscore <- LaplacesDemon::KLD(dens0$y, dens1$y)$sum.KLD.py.px
# KLscore
p.ls[[3]] <- ggplot(pparg, aes(x = icm_scores, fill = (surf_actives==1))) +
  geom_density(alpha = 0.5) + theme(legend.position="none") +
  annotate(geom="label", x=-350, y=0.04, label="ICM", hjust = 0) +
  annotate(geom="text", x=-350, y=0.03, label=paste("KL=",signif(KLscore,dig=3)), 
           hjust = 0, size=2) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6),
        axis.text=element_text(size=6))


# p-value plots -----------------------
p.ls[[4]] <- ~{
  par(mar=c(1.5,1.1,0.3,0.1), mgp=c(.5,.05,0), tck=-0.02, cex.lab = .5, cex.axis = .5)
  plot( 1, type="n", xlab="", ylab="p-value", xlim = c(1,ntested), 
        xaxt = 'n', ylim=c(1e-5, 1), log="y", cex.lab = .5, cex.axis = .5)
  mtext(side=1, line = .6, text = "Fraction Tested", cex = .5)
  ## axis(1, at = 1:ntested, labels = signif(tested/n,digits=1), cex.axis = .75)
  axis(1, at = c(1,3,6,8,10,13,15), labels=FALSE )
  ## labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"), cex.axis = .5)
  ## axis(3, at = 1:ntested, labels = tested, cex.axis = 0.4, tcl=NA)
  text(x = c(1,3,6,8,10,13,15),
       ## Move labels to just below bottom of chart.
       y =.3e-05, # par("usr")[3] ,
       ## Use names from the data list.
       labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       ## Adjust the labels to almost 100% right-justified.
       adj = 0.965,
       ## Increase label size.
       cex = .5)
  abline(h=.05,col="grey")
  for(k in 1:nmthd){
    lines( 1:ntested, pvals[1,,2,k], col=col.cd[k], lwd=1.5)
  }
}
p.ls[[5]] <- ~{
  par(mar=c(1.5,1.1,0.3,0.1), mgp=c(.5,.05,0), tck=-0.02, cex.lab = .5, cex.axis = .5)
  plot( 1, type="n", xlab="", ylab="p-value", xlim = c(1,ntested), 
        xaxt = 'n', ylim=c(1e-5, 1), log="y", cex.lab = .5, cex.axis = .5)
  mtext(side=1, line = .6, text = "Fraction Tested", cex = .5)
  ## axis(1, at = 1:ntested, labels = signif(tested/n,digits=1), cex.axis = .75)
  axis(1, at = c(1,3,6,8,10,13,15), labels=FALSE )
  ## labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"), cex.axis = .5)
  ## axis(3, at = 1:ntested, labels = tested, cex.axis = 0.4, tcl=NA)
  text(x = c(1,3,6,8,10,13,15),
       ## Move labels to just below bottom of chart.
       y =.3e-05, # par("usr")[3] ,
       ## Use names from the data list.
       labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       ## Adjust the labels to almost 100% right-justified.
       adj = 0.965,
       ## Increase label size.
       cex = .5)
  abline(h=.05,col="grey")
  for(k in 1:nmthd){
    lines( 1:ntested, pvals[1,,3,k], col=col.cd[k], lwd=1.5)
  }
}
p.ls[[6]] <- ~{
  par(mar=c(1.5,1.1,0.3,0.1), mgp=c(.5,.05,0), tck=-0.02, cex.lab = .5, cex.axis = .5)
  plot( 1, type="n", xlab="", ylab="p-value", xlim = c(1,ntested), 
        xaxt = 'n', ylim=c(1e-5, 1), log="y", cex.lab = .5, cex.axis = .5)
  mtext(side=1, line = .6, text = "Fraction Tested", cex = .5)
  ## axis(1, at = 1:ntested, labels = signif(tested/n,digits=1), cex.axis = .75)
  axis(1, at = c(1,3,6,8,10,13,15), labels=FALSE )
  ## labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"), cex.axis = .5)
  ## axis(3, at = 1:ntested, labels = tested, cex.axis = 0.4, tcl=NA)
  text(x = c(1,3,6,8,10,13,15),
       ## Move labels to just below bottom of chart.
       y =.3e-05, # par("usr")[3] ,
       ## Use names from the data list.
       labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       ## Adjust the labels to almost 100% right-justified.
       adj = 0.965,
       ## Increase label size.
       cex = .5)
  abline(h=.05,col="grey")
  for(k in 1:nmthd){
    lines( 1:ntested, pvals[2,,3,k], col=col.cd[k], lwd=1.5)
  }
}

# Plots of pointwise confidence intervals of difference of hit enrichment curves -----------------
p.ls[[7]] <- ~{
  par(mar=c(1.5,1.1,0.3,0.1), mgp=c(.5,.05,0), tck=-0.02, cex.lab = .5, cex.axis = .5)
  plot(1, type="n", xlab="", ylab="Difference", xaxt = 'n',
       xlim = c(1,ntested), ylim=c(-.2,.5), cex.lab = .5, cex.axis = .5)
  mtext(side=1, line = .6, text = "Fraction Tested", cex = .5)
  axis(1, at = c(1,3,6,8,10,13,15), labels=FALSE )
  ## labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"), cex.axis = .5)
  text(x = c(1,3,6,8,10,13,15), y = par("usr")[3]-.05,
       labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"),
       xpd = NA, srt = 35, adj = 0.965, cex = .5)
  abline(h=0,col="grey",lwd=2)
  for(k in 1:nmthd){
    lines( 1:ntested, recdiff[1,,2,k], col="black", lwd=1.5, lty=1 )
    lines( 1:ntested, cilow[1,,2,k], col=col.cd[k], lwd=1.5, lty=2 )
    lines( 1:ntested, cihigh[1,,2,k], col=col.cd[k], lwd=1.5, lty=2 )
  }
}
p.ls[[8]] <- ~{
  par(mar=c(1.5,1.1,0.3,0.1), mgp=c(.5,.05,0), tck=-0.02, cex.lab = .5, cex.axis = .5)
  plot(1, type="n", xlab="", ylab="Difference", xaxt = 'n',
       xlim = c(1,ntested), ylim=c(-.2,.5), cex.lab = .5, cex.axis = .5)
  mtext(side=1, line = .6, text = "Fraction Tested", cex = .5)
  axis(1, at = c(1,3,6,8,10,13,15), labels=FALSE )
  ## labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"), cex.axis = .5)
  text(x = c(1,3,6,8,10,13,15), y = par("usr")[3]-.05,
       labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"),
       xpd = NA, srt = 35, adj = 0.965, cex = .5)
  abline(h=0,col="grey",lwd=2)
  for(k in 1:nmthd){
    lines( 1:ntested, recdiff[1,,3,k], col="black", lwd=1.5, lty=1 )
    lines( 1:ntested, cilow[1,,3,k], col=col.cd[k], lwd=1.5, lty=2 )
    lines( 1:ntested, cihigh[1,,3,k], col=col.cd[k], lwd=1.5, lty=2 )
  }
}
p.ls[[9]] <- ~{
  par(mar=c(1.5,1.1,0.3,0.1), mgp=c(.5,.05,0), tck=-0.02, cex.lab = .5, cex.axis = .5)
  plot(1, type="n", xlab="", ylab="Difference", xaxt = 'n',
       xlim = c(1,ntested), ylim=c(-.2,.5), cex.lab = .5, cex.axis = .5)
  mtext(side=1, line = .6, text = "Fraction Tested", cex = .5)
  axis(1, at = c(1,3,6,8,10,13,15), labels=FALSE )
  ## labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"), cex.axis = .5)
  text(x = c(1,3,6,8,10,13,15), y = par("usr")[3]-.05,
       labels = c("0.001","0.002","0.005","0.01","0.02","0.05","0.1"),
       xpd = NA, srt = 35, adj = 0.965, cex = .5)
  abline(h=0,col="grey",lwd=2)
  for(k in 1:nmthd){
    lines( 1:ntested, recdiff[2,,3,k], col="black", lwd=1.5, lty=1 )
    lines( 1:ntested, cilow[2,,3,k], col=col.cd[k], lwd=1.5, lty=2 )
    lines( 1:ntested, cihigh[2,,3,k], col=col.cd[k], lwd=1.5, lty=2 )
  }
}

# Put all plots and legends together -----------------

# create plot grid
prow <- plot_grid(p.ls[[1]], p.ls[[4]], p.ls[[5]], 
                  p.ls[[7]], p.ls[[2]], p.ls[[6]], 
                  p.ls[[8]], p.ls[[9]], p.ls[[3]])

# add legends
df2 <- pparg
df2$surf_actives_labels <- ifelse(df2$surf_actives == 1, "+", "-")
p <- ggplot(df2, aes(x = maxz_scores, fill = surf_actives_labels)) +
  geom_density(alpha = 0.5) + 
  theme(legend.position = "bottom", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=14)) + 
  guides(fill=guide_legend(title="Activity Status:"))

legend_a <- get_legend(p)

df_tmp <- data.frame(pvals[i,,j,])
mthd <- c("EmProc", "CorrBinom", "McNemar", "IndJZ")
col.cd <- c("blue", "red", "cyan", "orange", "darkgreen")
colnames(df_tmp) <- mthd
df_tmp$x <- as.numeric(rownames(df_tmp))

df_tmp2 <- df_tmp %>% 
  pivot_longer(!x, names_to = "Method", values_to = "pval")

pp <- ggplot(df_tmp2, aes(x=x,y=pval, group=Method, color=Method)) +
  geom_line() +
  scale_colour_manual(values = col.cd[1:4],breaks=mthd[1:4]) + 
  labs(color="Procedures:") + 
  theme(legend.position = "bottom", 
        legend.title=element_text(size=10),
        legend.text=element_text(size=10)) 

legend_b <- get_legend(pp)

#legends <- plot_grid(legend_a, legend_b, nrow = 1)
legends <- plot_grid(legend_a, legend_b, nrow = 2)


final_plot <- plot_grid(prow, legends, ncol = 1, rel_heights = c( .95, .125))

# You may need to run this line twice to get the plot to show up in rstudio
# dont ask me why
pdf("Fig4.pdf",height = 6, width = 7)
final_plot
dev.off()