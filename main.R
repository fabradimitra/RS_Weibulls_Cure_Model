require(doParallel)
require(parallel)
require(nnet)
source("functions/sim.R")
source("functions/dgps.R")
source("functions/em.R")
source("functions/optimise_model_optim.R")
#
# Common parameters setting ----
tau <- 30 # 30 years follow-up: the follow up time in the data is in days!
lam0 <- 0.0367 # Scale baseline/expected surv
beta0 <- 5 # Shape baseline/expected surv
nsim <- 250  # number of simulation
maxit <- 10 # number of em iterations
# Simulation study ----
for (n in c(4000,8000,12000)){ 
  avcores <- parallel::detectCores()-1    # We always keep 1 spare core
  ncores <- min(n, avcores)               # We do not parallelize more than we need
  res <- mclapply(1:nsim, sim, 
                  n,tau,lam0,beta0,
                  components = 2,
                  maxit=10,
                  mc.cores = ncores)
  assign(paste("res",n,sep = ""), matrix(unlist(res), nrow = nsim, ncol = 3, byrow=TRUE),
         env =.GlobalEnv)
  #
  cat("End iters for n: ", n,"\n", file = "log/log.txt", append = TRUE)
}
#
save.image(file="data/data_mixwei2.RData") 
