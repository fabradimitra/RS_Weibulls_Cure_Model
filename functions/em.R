em <- function(t_obs, deltas, X, lam0, beta0,
               betad_start, gammas_start,
               pies_start,k=1,tol = 1e-6,maxit=1000,j){
  #
  n <- length(t_obs)
  pp1 <- dim(X)[2] # p + 1  ; covariates + column of ones
  #
  h0 <- beta0*lam0*(lam0*t_obs)^(beta0-1) # baseline hazard
  S0 <- exp(-(lam0*t_obs)^(beta0)) # baseline survival
  #
  log_lik <- function(deltas, probs, S0, h0, Sd, hd){
    fab <- (1-rowSums(probs))*S0*(h0^deltas)
    for(i in 1:k){
      fab <- fab + (probs[,i])*S0*Sd[,i]*((h0 + hd[,i])^(deltas))
    }
    fab[fab<1e-32]<-1e-32
    sum(log(fab))
  }
  #
  complete_loglik <- function(params, lb){
    k <- dim(params$Z)[2]
    ls <- matrix(data=0, nrow = params$n, ncol=k)
    out <- 0
    for (i in 1:k) {
      ls[,i] <- exp((params$X)%*%as.matrix(lb[((i-1)*params$pp1+1):(i*params$pp1)]))
      fab <- exp(-((ls[,i]*params$t_obs)^lb[k*params$pp1+i]))*
        ((params$h0+lb[k*params$pp1+i]*ls[,i]*
            ((ls[,i]*params$t_obs)^(lb[k*params$pp1+i]-1)))^params$deltas)
      fab[fab<1e-32]<-1e-32
      out <- out + sum(params$Z[,i]*log(fab))
    }
    return(out)
  }
  #
  Z <- matrix(data=0, nrow = n, ncol=k)
  pies <- pies_start
  gammas <- gammas_start
  betad <- betad_start
  betads <- matrix(betad_start,nrow = n, ncol = k, byrow = T)
  probs <- (exp(X%*%pies))/(1+rowSums(exp(X%*%pies)))
  lamds <- exp(X%*%gammas)
  #
  t_mat <- matrix(t_obs,nrow = n, ncol = k, byrow = F)
  deltas_mat <- matrix(deltas,nrow = n, ncol = k, byrow = F)
  h0s <-  matrix(h0,nrow = n, ncol = k, byrow = F)
  S0s <-  matrix(S0,nrow = n, ncol = k, byrow = F)
  #
  hds <- betads*lamds*((lamds*t_mat)^(betads-1)) # Excess hazard
  Sds <- exp(-(lamds*t_mat)^betads) # Excess survival
  # labels for the multinomial logit
  y<-c()
  for(i in 0:(k)){
    y<-c(y,rep(i,n))
  }
  #
  lliks<-c()
  log_lik_old <- log_lik(deltas, probs, S0, h0,
                         Sds, hds)
  iter <- 0
  dif <- 1
  cat("iter:", 0," dif:",dif," llk:", log_lik_old,"\n")
  #
  while ((dif > tol) & (iter<maxit)) {
    # Expectation
    w <- c()
    fab0 <- (1-rowSums(probs))*S0*(h0^deltas)
    den <- fab0
    for(i in 1:k){
      cac <- S0*Sds[,i]
      cac[cac<1e-16] <- 1e-16
      Z[,i] <- (probs[,i])*cac*((h0+hds[,i])^deltas)
      den <- den + Z[,i]
    }
    dens <- matrix(den, nrow = n, ncol = k, byrow = F)
    z0 <- fab0 / den
    Z <- Z/dens
    #
    # Maximization
    # Update pies
    w <- c(z0,as.vector(Z))
    w[w<1e-16] <- 1e-16
    X_glm <- X[,2:pp1]
    for(i in 1:k){
      X_glm <- rbind(X_glm,X[,2:pp1])
    }
    fab <- as.data.frame(cbind( X_glm, y))
    fab$y <- as.factor(fab$y)
    fab$y<-relevel(fab$y, ref="0") #ensure 0 is the reference category
    multinom_fit <- multinom(y~., data = fab, weights = w, trace=FALSE) # link logit
    if(k == 1){
      pies <- matrix(coef(multinom_fit),ncol = 1, nrow = pp1)
    }else{
      pies <- as.matrix(t(coef(multinom_fit)))
    }
    pies[which(is.na(pies))] <- 0
    probs <- (exp(X%*%pies))/(1+rowSums(exp(X%*%pies))) 
    prob0 <- 1-rowSums(probs)# probability of recovery
    #
    gammas_opt <- as.vector(gammas)
    res <- optim(par = c(gammas_opt,betad), 
                 fn = complete_loglik,
                 method = "Nelder-Mead",
                 params = list(t_obs=t_obs, Z=Z, deltas=deltas,
                               S0=S0, h0=h0, pp1 = pp1, X=X, n=n),
                 control=list(fnscale=-1))
    
    g_vec <- res$par[1:(k*pp1)]
    gammas <- matrix(g_vec, nrow = pp1, ncol = k, byrow = F)
    betad <- res$par[(k*pp1+1):(k*pp1+k)]
    betads <- matrix(betad,nrow = n, ncol = k, byrow = T)
    lamds <- exp(X%*%gammas)
    #
    hds <- betads*lamds*((lamds*t_mat)^(betads-1)) # Excess hazard
    Sds <- exp(-(lamds*t_mat)^betads) # Excess survival
    log_lik_cur <- log_lik(deltas, probs, S0, h0,
                           Sds, hds)
    # prepare for next iteration
    iter <- iter + 1
    lliks <- c(lliks,log_lik_cur)
    dif <- abs(log_lik_cur-log_lik_old)
    cat(j,"- components ",k, " - iter: ", iter," dif: ", dif," llk: ",log_lik_cur,"\n")
    log_lik_old <- log_lik_cur
  }
  return(list(lliks=lliks, pies = pies, gammas = gammas, 
              betad=betad, probs=probs, Z = Z, iter=iter, S0 = S0, h0 = h0))
} 