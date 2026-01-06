optimise_model_optim <- function(t_obs, deltas, X, 
                             S0, h0, betad_start,gammas_start, pies_start){
  #
  n <- length(t_obs)
  pp1 <- dim(X)[2] # p + 1
  k <- length(betad_start)
  pies <- as.vector(pies_start)
  gammas <- as.vector(gammas_start)
  betad <- betad_start
  #
  t_mat <- matrix(t_obs,nrow = n, ncol = k, byrow = F)
  #
  log_lik <- function(lb, params){
    #
    pies_mat <- matrix(lb[1:(params$k*params$pp1)], 
                        nrow = pp1, ncol = k, byrow = F)
    probs <- (exp(params$X%*%pies_mat))/(1+rowSums(exp(params$X%*%pies_mat)))
    #
    gammas_mat <- matrix(lb[(params$k*params$pp1+1):(2*params$k*params$pp1)], 
                          nrow = pp1, ncol = k, byrow = F)
    lamds <- exp(params$X%*%gammas_mat)
    #
    betad <- lb[(2*params$k*params$pp1+1):(2*params$k*params$pp1+params$k)]
    betad_mat <- matrix(betad, nrow = params$n, ncol = params$k, byrow = T)
    #
    hds <- betad_mat*lamds*((lamds*params$t_mat)^(betad_mat-1)) # Excess hazard
    Sds <- exp(-(lamds*params$t_mat)^betad_mat) # Excess survival
    # log-lik calculation
    fab <- (1-rowSums(probs))*params$S0*(params$h0^params$deltas)
    for(i in 1:params$k){
      fab <- fab + (probs[,i])*params$S0*Sds[,i]*((params$h0 + hds[,i])^(params$deltas))
    }
    fab[fab<1e-32]<-1e-32
    sum(log(fab))
  }                    
  res <- optim(par = c(pies_start,gammas_start,betad_start), 
               fn = log_lik,
               method = "Nelder-Mead",
               params = list(t_mat=t_mat, deltas=deltas,
                             S0=S0, h0=h0, pp1 = pp1, X=X, k=k, n=n),
               control=list(fnscale=-1),hessian = FALSE) 
  return(res)
}