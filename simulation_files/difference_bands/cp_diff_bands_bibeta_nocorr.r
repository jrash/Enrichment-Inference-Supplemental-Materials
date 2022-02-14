# Install my packages
list.of.packages <- c("doParallel","doMC","kableExtra", "foreach", "MASS", "copula", "foreach", "devtools", "KernSmooth")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

cl <- makeCluster(2)
registerDoParallel(cl)

setwd("Enrichment-Inference-Supplemental-Materials/simulation_files/difference_bands")
source("power_type1_cb.r")
source("hypothesis_test_new.r")

r.vec <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
r.vec <- r.vec[order(r.vec)]/150000
simreps = 10

method.v <- c("sup-t", "bonf", "EmProc")

dist <- "bibeta"
params.ls <- list(list(list(c(5, 2), c(2, 5)), list(c(4, 2), c(2, 5))),
                  list(list(c(5, 2), c(2, 5)), list(c(5, 2), c(2, 5))),
                  list(list(c(4, 2), c(2, 5)), list(c(4, 2), c(2, 5))))
rho <- c(.1, .1, .1)

k.ls <- rep(list(list()), length(params.ls)) 
k.c.ls <- rep(list(list()), length(params.ls)) 
for(i in seq_along(params.ls)) {
  for(j in seq_along(method.v)) {
    print(paste(params.ls[i], method.v[j]))
    k.ls[[i]][[j]] <- power.test(params= params.ls[[i]], dist = dist, pi.0.true = 1/500,
                                 myseed=569, method = method.v[j], plus = F,
                                 simreps = simreps, m = 150000,
                                 r.vec = r.vec, metric = "k", rho = rho[i])
    k.c.ls[[i]][[j]] <- power.test(params= params.ls[[i]], dist = dist, pi.0.true = 1/500,
                                   myseed=569, method = method.v[j], plus = T,
                                   simreps = simreps, m = 150000,
                                   r.vec = r.vec, metric = "k", rho = rho[i])
    save.image("power_band_bibeta_nocorr.rdata")
  }
}


