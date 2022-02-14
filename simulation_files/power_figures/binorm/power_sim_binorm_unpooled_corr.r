# Install my packages
list.of.packages <- c("doParallel","doMC","kableExtra", "foreach", "MASS", "copula", "foreach", "devtools", "KernSmooth")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# install_github("jrash/chemmodlab", force = T, quiet = T, upgrade = "never")
# library(chemmodlab)

cl <- makeCluster(6)
registerDoParallel(cl)

setwd("Enrichment-Inference-Supplemental-Materials/simulation_files/power_figures/binorm/")
source("../../power_type1.r")
source("../../hypothesis_test.r")

r.vec <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
r.vec <- r.vec[order(r.vec)]/150000
simreps = 10000

method.v <- c("AH", "binomial", "JZ ind", "mcnemar", "binomial ind")

dist <- "binorm"
params.ls <- list(list(c(sqrt(2) * .8, 0), c(sqrt(2) * .6, 0)), list(c(sqrt(2) * .8, 0), c(sqrt(2) * .8, 0)),
                  list(c(sqrt(2) * .6, 0), c(sqrt(2) * .6, 0)))
rho <- c(.9, .9, .9)

k.ls <- rep(list(list()), length(params.ls)) 
k.c.ls <- rep(list(list()), length(params.ls)) 
for(i in seq_along(params.ls)) {
  for(j in seq_along(method.v)) {
    print(paste(params.ls[i], method.v[j]))
    k.ls[[i]][[j]] <- power.test(params= params.ls[[i]], dist = dist, pi.0.true = 1/500,
                                 myseed=569, method = method.v[j], plus = F, pool = F, 
                                 simreps = simreps, m = 150000,
                                 r.vec = r.vec, metric = "k", rho = rho[i])
    k.c.ls[[i]][[j]] <- power.test(params= params.ls[[i]], dist = dist, pi.0.true = 1/500,
                                   myseed=569, method = method.v[j], plus = T, pool = F,
                                   simreps = simreps, m = 150000,
                                   r.vec = r.vec, metric = "k", rho = rho[i])
    save.image("power_all_binorm_unpooled_corr.rdata")
  }
}


