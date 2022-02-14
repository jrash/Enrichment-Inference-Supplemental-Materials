
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

# # # Simulate Case 1 scores from paper
# r.vec <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
# r.vec <- r.vec[order(r.vec)]/150000
# 
# m <- 150000
# pi.0.true <- 1/501
# set.seed(111)
# X <- rbinom(n=m, size=1, prob=pi.0.true)
# mu.mix <- c(1.4, 0)
# sd.mix <- 1
# S <- rnorm(n=m, mu.mix[1], sd.mix)*X + rnorm(n=m, mu.mix[2], sd.mix)*(1-X)
# 
# metric = "k"; type = "band"; method = "sup-t"; correction = "JZ";
# conf.level = .95; boot.rep = 100; myseed = 111; mc.rep = 100000;
# simreps = 10; pi.0.true = 1/501;dist = "binorm"; plus2 = T
# 
# params <- c(1.4, 0)


# params= params.ls[[i]]; dist = dist.v[i]; pi.0.true = 1/501;
# myseed=567; method = "JZ"; plus2 = correction.v[j];
# simreps = simreps; m = 150000; B = 1000; mc.rep = mc.rep;
# r.vec = r.vec; metric = "k"

# Compare to simulated interval
cov.prob.tab = function(params, pi.0.true, myseed=567, dist = "binorm",
                    simreps = 100, r.vec = c(.01, .1, .5), metric = "k",
                    m = 100, plus = F, method = "JZ", B = 100, mc.rep = 1000000){
  r <- (1:m)/m
  idx <- which(r %in% r.vec)
  r <- r[idx]
  set.seed(myseed)
  
  if(dist == "binorm" | dist == "bibeta") {
    suppressWarnings(c <- qmix(1-r, params, c(pi.0.true, 1-pi.0.true), expand=.1, dist = dist))
  } else if(dist == "unif") {
    p <- 1-r
    # mixture weights
    w1 <- (1-pi.0.true)
    w2 <- (pi.0.true)
    # knots
    t1 <- w1*(.25/.75)
    t2 <- w1 + w2*((.75-.25)/.75)
    # p <- runif(1000000)
    c <- vector(length = length(p))
    for(j in seq_along(p)) {
      if(p[j] <= 0) {
        c[j] <- 0
      } else if(p[j] > 0 & p[j] <= t1) {
        c[j] <- .75*p[j]/w1
      } else if(p[j] > t1 & p[j] <= t2) {
        c[j] <- .75*p[j] + w2*.25
      } else if(p[j] > t2 & p[j] <= 1) {
        c[j] <- .75*(p[j]-w1)/w2 + .25
      } else if(p[j] == 1) {
          c[j] <- 1L
      }
    }
  }

  if (dist == "binorm") {
    sd.s = 1
    F.plus <- pnorm(c, mean = params[1], sd.s)
  } else if (dist == "bibeta") {
    F.plus <- pbeta(c, params[[1]][1], params[[1]][2])
  } else if (dist == "unif") {
    F.plus <- punif(c, params[[1]][1], params[[1]][2])
  }
  
  k.true <- vector(length = length(F.plus))
  pi.true <- vector(length = length(F.plus))
  
  # if k or pi is above 1, then they should actually be slightly below 1
  k.true <- ifelse(k.true > 1, 1 - 1E-10, k.true)
  pi.true <- ifelse(pi.true > 1, 1 - 1E-10, pi.true)

  for (i in 1:length(k.true)) {
    # Account for the case when you have non-overlapping uniforms and you 
    # know that pi should be 1
    if (dist == "unif" & c[i] > params[[2]][2]) {
      pi.true[i] <- 1L
      k.true[i] <- r[i]/pi.0.true
    } else {
      k.true[i] <- 1 - F.plus[i]
      pi.true[i] <- k.true[i]*pi.0.true/r[i]
    }
  }
  eff_seeds <- sample(1:2^15, simreps)
  results <- foreach(z = 1:simreps, .export=c("mvrnorm", "PerfCurveBands", "EstLambda")) %dopar% {
    set.seed(eff_seeds[z])
    X = rbinom(n=m, size=1, prob=pi.0.true)
    if(dist == "binorm") {
      S = rnorm(n=m, params[1], sd.s)*X + rnorm(n=m, params[2], sd.s)*(1-X)
    } else if(dist == "bibeta") {
      S = rbeta(n=m, params[[1]][1], params[[1]][2])*X + rbeta(n=m, params[[2]][1], params[[2]][2])*(1-X)
    } else if(dist == "unif") {
      S = runif(n=m, params[[1]][1], params[[1]][2])*X + runif(n=m, params[[2]][1], params[[2]][2])*(1-X)
    } 
    
    PerfCurveBands(S, X, r, metric = "rec", type = ifelse(method == "JZ", "pointwise", "band"), method = method, plus = plus, conf.level = .95,
                   mc.rep = mc.rep, myseed = eff_seeds[z])
  }
  
  covout <- matrix(NA, nrow=1, ncol=simreps)
  widthout <- matrix(NA, nrow=length(r), ncol=simreps)
  interval.ls <- list()
  params <- matrix(NA, nrow=length(r), ncol=simreps)
  rec.mat <- matrix(NA, nrow=length(r), ncol=simreps)
  
  for (i in 1:simreps) {
    ind <- vector(length = length(k.true))
    for(j in seq_along(k.true)) {
      ind[j] <- (round(results[[i]]$CI[j, 1], 10) <= k.true[j]) & (round(results[[i]]$CI[j, 2], 10) >= k.true[j])
    }
    covout[i] <- all(ind)
    widthout[, i] <- results[[i]]$CI[, 2] - results[[i]]$CI[, 1] 
    interval.ls[[i]] <- results[[i]]$CI
    params[, i] <- k.true
    rec.mat[, i] <- results[[i]]$rec
  }
  
  return(list(r = r, covout = covout, widthout = widthout, interval.ls = interval.ls, params = params, rec.mat = rec.mat))
}

