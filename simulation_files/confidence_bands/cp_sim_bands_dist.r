# Install my packages
list.of.packages <- c("doParallel","kableExtra", "foreach", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

cl <- makeCluster(2)
registerDoParallel(cl)

setwd("Enrichment-Inference-Supplemental-Materials/simulation_files/confidence_bands/")
source("coverage_probabilities_bands.r")
source("confidence_bands.r")

dist.v <- c("binorm", "binorm",
            "bibeta", "bibeta",
            "unif")

r.vec <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
r.vec <- r.vec[order(r.vec)]/150000

params.ls <- list(c(1.4, 0), c(.5, 0), list(c(5,2), c(2,5)),
                  list(c(20,1), c(1,20)), list(c(.25,1), c(0,.75)))
method.v <- c("sup-t", "sup-t", "bonf", "bonf", "JZ", "JZ")
correction.v <- c(F, T, F, T, F, T)
simreps <- 5
mc.rep <- 100000

k.ls <- list()
pi.ls <- list()
p <- length(method.v)
i <- 3; j <- 2
for(i in seq_along(dist.v)) {
  for(j in seq_along(method.v)) {
    print(paste(method.v[j]))
    k.ls[[(i-1) * p + j]] <- cov.prob.tab(params= params.ls[[i]], dist = dist.v[i], pi.0.true = 1/501,
                                myseed=567, method = method.v[j], plus = correction.v[j],
                                simreps = simreps, m = 150000, B = 1000, mc.rep = mc.rep,
                                r.vec = r.vec, metric = "k")
  }
  save.image("cp_all_bands_dist.rdata")
}


