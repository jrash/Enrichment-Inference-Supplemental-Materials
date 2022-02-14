
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
  
  
  return(k.true1 - k.true2)
}


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
