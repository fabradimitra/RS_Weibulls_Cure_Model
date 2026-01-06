dgps<-function(n, X, tau, probs, lamd1, lamd2, 
               beta1, beta2, lam0, beta0, maxtry=100){
  
  # n: number of observations
  # X: matrix of covariates' values for the n observations and an initial column of ones
  # tau: maximum follow-up time in days
  # probs_vec: vector of probability of being cured for each individual
  # lamd_vec: vector of excess hazard of each individual
  # betad: Shape excess surv
  # lam0: Scale baseline/expected surv
  # beta0: Shape baseline/expected surv
  # cp: censoring proportion
  env_dgps<-new.env()
  # z=0 cured, z=1 uncured group 1, z=2 uncured group 2
  z <- apply(probs, MARGIN=1, FUN = function(p){
    sample(x=c(0,1,2), size = 1, prob = p, replace = TRUE)
  })
  n0 <- length(which(z==0))
  n1 <- length(which(z==1))
  n2 <- length(which(z==2))
  # Generate time to event for cured proportion 
  t_cure <- rweibull(n0, scale = (1/lam0), shape = beta0)
  # Generate time to event for uncured proportion 
  # t_uncure_1 <- rweibull(n1, scale = (1/lamd1), shape = beta1)
  # t_uncure_2 <- rweibull(n2, scale = (1/lamd2), shape = beta2)
  u <- runif(n1)
  finvuc <- function(x, u, lamd, betad){
    u-exp(-(lam0*x)^(beta0)-(lamd*x)^(betad))
  }
  t_uncure_1 <- rep(NA, n1)
  lamd_sel <- lamd1[z==1] # Select unit-specific hazard
  assign("go", FALSE, env=env_dgps)
  assign("ntry", 1, env=env_dgps)
  assign("x", 0, env=env_dgps)
  for (i in 1:n1){
    assign("unf", u[i], env=env_dgps)
    while ((!get("go", env=env_dgps))&(get("ntry", env=env_dgps)<maxtry)) {
      tryCatch(
        {
          assign("x",uniroot(finvuc, interval = c(0,tau),
                             extendInt = "upX", tol=1e-12,
                             u = get("unf", env=env_dgps),
                             lamd = lamd_sel[i], betad = beta1)$root,
                 env=env_dgps)
          assign("go", TRUE, env=env_dgps)
          
        },
        error = function(e){
          assign("unf", runif(1), env=env_dgps)
          assign("ntry", get("ntry", env=env_dgps)+1, env=env_dgps)
        }
      )
    }
    t_uncure_1[i] <- get("x", env=env_dgps)
    assign("go", FALSE, env=env_dgps)
    assign("ntry", 1, env=env_dgps)
  }
  t_uncure_1[t_uncure_1<1e-32] <- 1e-32
  u <- runif(n2)
  t_uncure_2 <- rep(NA, n2)
  lamd_sel <- lamd2[z==2] # Select unit-specific hazard
  assign("go", FALSE, env=env_dgps)
  assign("ntry", 1, env=env_dgps)
  assign("x", 0, env=env_dgps)
  for (i in 1:n2){
    assign("unf", u[i], env=env_dgps)
    while ((!get("go", env=env_dgps))&(get("ntry", env=env_dgps)<maxtry)) {
      tryCatch(
        {
          assign("x",uniroot(finvuc, interval = c(0,tau),
                             extendInt = "upX", tol=1e-12,
                             u = get("unf", env=env_dgps),
                             lamd = lamd_sel[i], betad = beta2)$root,
                 env=env_dgps)
          assign("go", TRUE, env=env_dgps)
          
        },
        error = function(e){
          assign("unf", runif(1), env=env_dgps)
          assign("ntry", get("ntry", env=env_dgps)+1, env=env_dgps)
        }
      )
    }
    t_uncure_2[i] <- get("x", env=env_dgps)
    assign("go", FALSE, env=env_dgps)
    assign("ntry", 1, env=env_dgps)
  }
  t_uncure_2[t_uncure_2<1e-32] <- 1e-32
  t_tot <- rep(NA,n)
  t_tot[z==0]<-t_cure
  t_tot[z==1]<-t_uncure_1
  t_tot[z==2]<-t_uncure_2
  # Administrative Censoring
  deltas <- rep(1,n)
  deltas[which(t_tot>tau)] <- 0
  t_tot[deltas==0] <- tau
  t_tot[t_tot<1e-32]<-1e-32
  return(list(t_tot=t_tot,deltas=deltas,z=z))
} 