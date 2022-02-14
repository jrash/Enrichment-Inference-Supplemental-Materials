
power.test = function(params, pi.0.true, myseed=567, dist = "binorm",
                      simreps = 100, r.vec = c(.01, .1, .5), metric = "k",
                      m = 100, plus = T, pool, method = "AH",
                      rho, s.vec.ls = list(c(1, 1), c(1, 1))) {
  
  eff_seeds <- sample(1:2^15, simreps)
  results <- foreach(z = 1:simreps, .export=c("normalCopula", "mvdc", "rMvdc", "PerfCurveTest", "mvrnorm", "EstLambda")) %dopar% {
    set.seed(eff_seeds[z])
    X = rbinom(n=m, size=1, prob=pi.0.true)
    if (sum(X)==0) X = rbinom(n=m, size=1, prob=pi.0.true)
    if(dist == "binorm") {
      # Target parameters for univariate normal distributions
      # Parameters for bivariate normal distribution
      s1 <- s.vec.ls[[1]][1] ; s2 <- s.vec.ls[[2]][1];
      sigma1 <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2)
      s1 <- s.vec.ls[[1]][2]; s2 <- s.vec.ls[[2]][2];
      sigma2 <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2)# Covariance matrix
      bvn1 <- mvrnorm(m, mu = c(params[[1]][1], params[[2]][1]), Sigma = sigma1)*X + mvrnorm(m, mu = c(params[[1]][2], params[[2]][2]), Sigma = sigma2)*(1-X)
      S1 <- bvn1[, 1]
      S2 <- bvn1[, 2]
    } else if(dist == "bibeta") {
      myCopPlus <- normalCopula(param=c(rho), dim = 2, dispstr = "un")
      myMvdPlus <- mvdc(copula=myCopPlus, margins=c("beta", "beta"),
                        paramMargins=list(list(shape1=params[[1]][[1]][1], shape2=params[[1]][[1]][2]),
                                          list(shape1=params[[2]][[1]][1], shape2=params[[2]][[1]][2])))
      myCopNeg <- normalCopula(param=c(rho), dim = 2, dispstr = "un")
      myMvdNeg <- mvdc(copula=myCopNeg, margins=c("beta", "beta"),
                       paramMargins=list(list(shape1=params[[1]][[2]][1], shape2=params[[1]][[2]][2]),
                                         list(shape1=params[[2]][[2]][1], shape2=params[[2]][[2]][2])))
      
      bvn1 <- rMvdc(m, myMvdPlus)*X + rMvdc(m, myMvdNeg)*(1-X)
      S1 <- bvn1[, 1]
      S2 <- bvn1[, 2]
    } else if(dist == "unif") {
      myCopPlus <- normalCopula(param=c(rho), dim = 2, dispstr = "un")
      myMvdPlus <- mvdc(copula=myCopPlus, margins=c("unif", "unif"),
                        paramMargins=list(list(min=params[[1]][[1]][1], max=params[[1]][[1]][2]),
                                          list(min=params[[2]][[1]][1], max=params[[2]][[1]][2])))
      myCopNeg <- normalCopula(param=c(rho), dim = 2, dispstr = "un")
      myMvdNeg <- mvdc(copula=myCopNeg, margins=c("unif", "unif"),
                       paramMargins=list(list(min=params[[1]][[2]][1], max=params[[1]][[2]][2]),
                                         list(min=params[[2]][[2]][1], max=params[[2]][[2]][2])))
      
      bvn1 <- rMvdc(m, myMvdPlus)*X + rMvdc(m, myMvdNeg)*(1-X)
      S1 <- bvn1[, 1]
      S2 <- bvn1[, 2]
    } 
    
    PerfCurveTest(S1, S2, X, r = r.vec, metric = "rec", method = method, plus = plus, pool = pool, alpha = .05, seed = eff_seeds[z])
  }
  
  sr <- length(results)
  covout <- matrix(NA, nrow=length(r.vec), ncol=sr)
  widthout <- matrix(NA, nrow=length(r.vec), ncol=sr)
  interval.ls <- list()
  
  for (i in 1:sr) {
    covout[, i] <- results[[i]]$p_value < .05
    widthout[, i] <- results[[i]]$ci_interval[, 2] - results[[i]]$ci_interval[, 1] 
    interval.ls[[i]] <- results[[i]]$ci_interval
  }
  
  return(list(r = r.vec, covout = covout, widthout = widthout, interval.ls = interval.ls, params = params, eff_seeds))
  
}