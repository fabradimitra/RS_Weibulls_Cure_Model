sim <- function(j,n,tau,lam0,beta0,components,tol=1e-8,maxit=100){
  tryCatch( 
    { 
      set.seed(j)
      # Generate some covariates
      ones <- rep(1, n)
      # Age
      x1 <-rpois(n, 65)
      #hist(x1)
      # Period of diagnosis
      x2_1 <- runif(n/3, 0, tau/2) 
      x2_2 <- runif(ceiling(2*n/3), tau/2, tau) # In later years more screenings
      x2 <- c(x2_1,x2_2)
      # Sex: 1 female, 0 male
      x3 <- rbinom(n, 1, 0.48) 
      # Put variables together
      X <- cbind(ones, scale(x1, center = FALSE), scale(x2, center = FALSE), x3)
      colnames(X) <- c("intercept","age","period_of_diagnosis", "sex")
      pp1 <- dim(X)[2]
      #
      # Generate true parameters
      pies1_true <- c(runif(1, min = -1.0, max = -0.6),    # Intercept (affects Component 1 probability)
                      runif(1, min = 0.2, max = 0.4),     # Age coefficient (affects Component 2)
                      runif(1, min = -0.6, max = -0.2),   # Period of diagnosis (affects Component 1)
                      runif(1, min = -0.2, max = 0.2))    # Sex coefficient (small effect)     
      #
      pies2_true <- c(runif(1, min = -1.5, max = -1.0),   # Intercept (affects Component 0 probability)
                      runif(1, min = -0.5, max = 0.2),    # Age coefficient (affects Component 1)
                      runif(1, min = -0.5, max = -0.3),   # Period of diagnosis (affects Component 2)
                      runif(1, min = -0.2, max = 0.2))    # Sex coefficient (small effect)
      # Generate Data
      den <- 1+exp(X%*%pies1_true)+ exp(X%*%pies2_true)
      probs1 <- (exp(X%*%pies1_true))/den
      probs2 <- (exp(X%*%pies2_true))/den
      probs0 <- 1-probs1-probs2
      probs <- cbind(probs0,probs1,probs2)
      #
      gammas1_true <- c(
        runif(1, min = -2.5, max = -2),  # Intercept: log(0.178) ≈ -1.73
        runif(1, min = 0.05, max = 0.2),   # Small age effect
        -runif(1, min = 0.1, max = 0.3),    # Small period effect
        runif(1, min = -0.2, max = 0.2)     # Negligible sex effect
      )
      
      # Weibull parameters for Component 2 (Mean ~10)
      gammas2_true <- c(
        runif(1, min = -0.8, max = -0.6),  # Intercept: log(0.481) ≈ -0.732
        runif(1, min = 0.1, max = 0.3),    # Small age effect
        -runif(1, min = 0.3, max = 0.5),    # Moderate period effect
        runif(1, min = -0.2, max = 0.2)     # Negligible sex effect
      )    
      #
      # Compute lambda values
      lamd1_true <- exp(X %*% gammas1_true)
      lamd2_true <- exp(X %*% gammas2_true)
      #
      beta1_true <- runif(1, min = 1.5, max = 2)  
      beta2_true <- runif(1, min = 1, max = 1.1)    
      #
      dat <-dgps(n, X, tau, probs, lamd1_true, lamd2_true, 
                 beta1_true, beta2_true, lam0, beta0)
      #
      if(components==-1){
        components <- 1
        bic_prev <- Inf
        bic_cur <- Inf
        while(bic_cur<=bic_prev){
          bic_prev <- bic_cur
          res_em <-em(t_obs = dat$t_tot,
                      deltas = dat$deltas,
                      X = X,
                      lam0 = lam0,
                      beta0 = beta0,
                      betad_start = rep(1,components),
                      gammas_start = matrix(rep(1, dim(X)[2]),
                                            nrow = pp1, ncol = components,
                                            byrow = F),
                      pies_start = matrix(rep(1, dim(X)[2]),
                                          nrow = pp1, ncol = components,
                                          byrow = F),
                      k=components,
                      tol = 1e-8,
                      maxit=maxit,
                      j=j)
          #
          res <- optimise_model_optim(t_obs = dat$t_tot,
                                      deltas = dat$deltas,
                                      X = X,
                                      S0 = res_em$S0, 
                                      h0 = res_em$h0,
                                      betad_start = res_em$betad,
                                      gammas_start = res_em$gammas,
                                      pies_start = res_em$pies)
          bic_cur <- length(res$par)*n - res$value
          components <- components + 1
          cat(j,"- components ",components, " - parameters estimated" ,"\n")
        }
        components <- components - 1
      }
      res_em <-em(t_obs = dat$t_tot,
                  deltas = dat$deltas,
                  X = X,
                  lam0 = lam0,
                  beta0 = beta0,
                  betad_start = rep(1,components),
                  gammas_start = matrix(rep(1, dim(X)[2]),
                                        nrow = pp1, ncol = components,
                                        byrow = F),
                  pies_start = matrix(rep(1, dim(X)[2]),
                                      nrow = pp1, ncol = components,
                                      byrow = F),
                  k=components,
                  tol = 1e-8,
                  maxit=maxit,
                  j=j)
      #
      res <- optimise_model_optim(t_obs = dat$t_tot,
                                  deltas = dat$deltas,
                                  X = X,
                                  S0 = res_em$S0, 
                                  h0 = res_em$h0,
                                  betad_start = res_em$betad,
                                  gammas_start = res_em$gammas,
                                  pies_start = res_em$pies)
      pies <- matrix(res$par[1:(components*pp1)], 
                     nrow = pp1, ncol = components, byrow = F)
      #
      pies[which(is.na(pies))] <- 0
      probs <- (exp(X%*%pies))/(1+rowSums(exp(X%*%pies))) 
      prob0_est <- 1-rowSums(probs) # estimated cure proportions
      cat(j, ": ", mean((probs0 - prob0_est)^2), ", ", 
          mean(abs(probs0 - prob0_est)),", ", components, "\n",
          file="log/log.txt", append = TRUE)
      #
      write("\n",file="log/log.txt",append=TRUE)
      #
      return(list(prob0_MAE = mean(abs(probs0 - prob0_est)),
                  prob0_MSE = mean((probs0 - prob0_est)^2),
                  n_comp = components)
      )
    },
    error = function(w){
      cat("Error at iteration:",j, 
          file="log/log.txt", sep = " ", append = TRUE)
      write("\n",file="log/log.txt",append=TRUE)
      cat("Error:",warning(w), 
          file="log/log.txt", sep = " ", append = TRUE)
      write("\n",file="log/log.txt",append=TRUE)
    }
  )
  #
}
