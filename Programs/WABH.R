### This function computes the FDR constraint in *
FDR.constraint<-function(lambda, g, p,alpha = 0.05){
  #Calculate the threshold on the log scale
  c_m <- lambda*(1-p)*(1-alpha)/(p*(1+lambda*alpha))
  log_t<- pnorm(log(c_m)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  log_t[log_t>0] <- 0
  t <- exp(log_t)
  #Calculate the log-power then the power.
  log_pow <- pnorm(qnorm(log_t,lower.tail = FALSE, log.p = TRUE)-g, lower.tail =FALSE, log.p = TRUE)
  pow <- exp(log_pow)
  #Calculate the expected FDR.
  efdr <- mean(t*(1-p))/(mean(t*(1-p)) + mean(p*pow))
  return((log(efdr) - log(alpha))^2)
}

### This uses optimize to find the value of lambda
### that satisfies the constraint above.  May have to mess with
### the endpoints.  Theoretically in (0, infinity).
get.lambda<-function(g, p, alpha=.05){
  FDR.c<-function(lambda)
  {
    FDR.constraint(lambda,g,p, alpha)
  }
  
  max_v <- 10
  t1 <- try(asd <- optimize( FDR.c, c(1e-100, max_v), maximum = FALSE), silent = TRUE)
  lambda <- asd$minimum
  while(lambda > max_v*0.95 & max_v < 1e100){
    max_v <- max_v*10
    t1 <- try(asd <- optimize( FDR.c, c(1e-100, max_v), maximum = FALSE), silent = TRUE)
    lambda <- asd$minimum
  }
  
  if(lambda > max_v*0.95 & max_v >= 1e100){stop(t1[1],"Error: value lambda reach 1e100.")}
  return(lambda)
}

### This function computes the FDR constraint in *
FDR.constraint.eta<-function(lambda, s, p,alpha = 0.05,p_MMW){
  #This function doesn't require eta, so first eta and "g" are calculated.
  c_MMW <- lambda*(1-p_MMW)*(1-alpha)/(p_MMW*(1+lambda*alpha))
  eta  <- sqrt(2*log(max(c(c_MMW,exp(1e-5)))))*min(s)
  g    <- eta/s
  
  #Calculate the threshold on the log scale
  c_m <- lambda*(1-p)*(1-alpha)/(p*(1+lambda*alpha))
  log_t<- pnorm(log(c_m)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  log_t[log_t>0] <- 0
  t <- exp(log_t)
  #Calculate the log-power then the power.
  log_pow <- pnorm(qnorm(log_t,lower.tail = FALSE, log.p = TRUE)-g, lower.tail =FALSE, log.p = TRUE)
  pow <- exp(log_pow)
  #Calculate the expected FDR.
  efdr <- mean(t*(1-p))/(mean(t*(1-p)) + mean(p*pow))
  return((log(efdr) - log(alpha))^2)
}


### This uses optimize to find the value of lambda
### that satisfies the constraint above.  May have to mess with
### the endpoints.  Theoretically in (0, infinity).
get.lambda.eta<-function(s, p, alpha=.05,p_MMW){
  FDR.c<-function(lambda)
  {
    FDR.constraint.eta(lambda,s,p, alpha,p_MMW)
  }
  
  max_v <- 10
  t1 <- try(asd <- optimize( FDR.c, c(1e-5, max_v), maximum = FALSE), silent = TRUE)
  lambda <- asd$minimum
  while(lambda > max_v*0.95 & max_v < 1e100){
    max_v <- max_v*10
    t1 <- try(asd <- optimize( FDR.c, c(1e-5, max_v), maximum = FALSE), silent = TRUE)
    lambda <- asd$minimum
  }
  
  if(lambda > max_v*0.95 & max_v >= 1e100){stop(t1[1],"Error: MMW criteria failed. Consider changing p_MMW.\n")}
  return(lambda)
}

## This puts it all together to return weights the most important function in the paper

get.logtm <- function(g,p,alpha){
  lambda<-get.lambda(g,p,alpha)
  c_m <- lambda*(1-p)*(1-alpha)/(p*(1+lambda*alpha))
  log_t<- pnorm(log(c_m)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  log_t[log_t>0] <- 0
  t <- exp(log_t)
  return(list(log_t=log_t,lambda=lambda))
}

get.logtm.eta <- function(s,p,alpha,p_MMW){
  lambda <-get.lambda.eta(s,p,alpha,p_MMW)
  c_MMW <- lambda*(1-p_MMW)*(1-alpha)/(p_MMW*(1+lambda*alpha))
  if(c_MMW<exp(1e-5)){cat("Warning: MMW criteria failed. Likely due to there being too much power to have monotone weights. \n  Consider changing p_MMW.\n")}
  eta  <- sqrt(2*log(max(c(c_MMW,exp(1e-5)))))*min(s)
  
  g    <- eta/s
  
  #Calculate the threshold on the log scale
  c_m <- lambda*(1-p)*(1-alpha)/(p*(1+lambda*alpha))
  log_t<- pnorm(log(c_m)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  log_t[log_t>0] <- 0
  t <- exp(log_t)
  return(list(log_t=log_t,lambda=lambda,eta=eta))
}


get.weights <-function(g, p, alpha,p_MMW,find_eta=FALSE){
  if(!find_eta){
    t_find <- get.logtm(g,p,alpha)
    thresholds <- exp(t_find$log_t)
    eta=NULL
    lambda = t_find$lambda
  }
  if(find_eta){
    t_eta_find <- get.logtm.eta(g,p,alpha,p_MMW)
    thresholds <- exp(t_eta_find$log_t)
    
    eta=t_eta_find$eta   
    lambda = t_eta_find$lambda 
  }
  return(list(weights = thresholds/mean(thresholds),lambda=lambda,eta=eta))
}




weighted_p <- function(p_vals,pred_SE, pi_hat, alpha, eta, p_MMW){
  #Get two versions of the weighted p-values from pi_hat, eta.
  
  ### Calculating the weights when using the MMW criteria for eta.
  weights_info  <- get.weights(g = pred_SE, p = pi_hat, alpha=alpha, p_MMW, find_eta = TRUE)
  
  ### Calculating the weights when eta is given.
  eta_pow <- eta #*0.4
  weights2 <- get.weights(g = eta_pow/pred_SE, p = pi_hat, alpha = alpha, p_MMW, find_eta = FALSE)$weights
  weights <- weights_info$weights
  
  
  p_weight <- p_vals/weights
  p_weight[p_weight>1] <- 1
  p_weight2 <- p_vals/weights2
  p_weight2[p_weight2>1] <- 1
  
  return(list(p_weight_MMW=p_weight, p_weight_eta = p_weight2, eta_MMW = weights_info$eta, MMW_weights = weights_info$weights, eta_weights =weights2 ))
  
}

MTR <- function(p_vals,p_wbonf,p_wbh,X_mean2,alpha,p_hat=1,signal=NULL,per_val=0.1){
  
  alpha_hat <- alpha/p_hat
  p_vals2 <- p_vals
  p_vals2[X_mean2<=per_val | X_mean2 >= (1-per_val)] <- 1
  
  W_BH <- p.adjust(p_wbonf,method = "BH")
  W_BH2<- p.adjust(p_wbh,method = "BH")
  T_R <- p_vals2
  R2_BH <- p.adjust(p_vals2[X_mean2>per_val & X_mean2 < (1-per_val)],method = "BH")
  T_R[X_mean2>per_val & X_mean2 < (1-per_val)] <- R2_BH
  R_BH <- p.adjust(p_vals,method = "BH")
  
  BH_corr_pvals <- data.frame(WABH = W_BH, WABH_2 = W_BH2, Ten = T_R,BH = R_BH,ABH = R_BH)
  
  BH_res <- data.frame(WABH = 1*I(W_BH<alpha_hat), WABH_2 = 1*I(W_BH2<alpha_hat), Ten = 1*I(T_R<alpha),BH = 1*I(R_BH<alpha),ABH = 1*I(R_BH<alpha_hat))
  
  W_BY <- p.adjust(p_wbonf,method = "holm")
  W_BY2<- p.adjust(p_wbh,method = "holm")
  R_BY <- p.adjust(p_vals,method = "holm")
  T_R_Bonf <- p_vals2
  R2_BY <- p.adjust(p_vals2[X_mean2>per_val & X_mean2 < (1-per_val)],method = "holm")
  T_R_Bonf[X_mean2>per_val & X_mean2 < (1-per_val)] <- R2_BY
  
  Bonf_res <- data.frame(WABH = 1*I(W_BY<alpha), WABH_2 = 1*I(W_BY2<alpha), Ten = 1*I(T_R_Bonf<alpha),Bonf = 1*I(R_BY<alpha))
  
  BH_sum <- NULL
  Bonf_sum <- NULL
  if(!is.null(signal)){
    M_new <- length(X_mean2)
    n_signal <- !(1:M_new %in% signal)
    BH_sum   <- c(sum(BH_res$BH),sum(BH_res$BH[signal]),sum(BH_res$BH[n_signal]))
    W_res    <- c(sum(BH_res$WABH),sum(BH_res$WABH[signal]),sum(BH_res$WABH[n_signal]))
    W2_res   <- c(sum(BH_res$WABH_2),sum(BH_res$WABH_2[signal]),sum(BH_res$WABH_2[n_signal]))
    ADAP_res <- c(sum(BH_res$ABH),sum(BH_res$ABH[signal]),sum(BH_res$ABH[n_signal]))
    Ten_res  <- c(sum(BH_res$Ten),sum(BH_res$Ten[signal]),sum(BH_res$Ten[n_signal]))
    
    BH_sum <- data.frame(WABH_sum = W_res, WABH_sum_2 = W2_res, Ten_sum = Ten_res, ABH_sum = ADAP_res,BH_sum = BH_sum)
    
    W_res2   <- c(sum(Bonf_res$WABH),sum(Bonf_res$WABH[signal]),sum(Bonf_res$WABH[n_signal]))
    W2_res2  <- c(sum(Bonf_res$WABH_2),sum(Bonf_res$WABH_2[signal]),sum(Bonf_res$WABH_2[n_signal]))
    Bonf_sum <- c(sum(Bonf_res$Bonf),sum(Bonf_res$Bonf[signal]),sum(Bonf_res$Bonf[n_signal]))
    Ten_res2 <- c(sum(Bonf_res$Ten),sum(Bonf_res$Ten[signal]),sum(Bonf_res$Ten[n_signal]))
    
    Bonf_sum <- data.frame(WABH_sum = W_res2, WABH_sum_2 = W2_res2, Ten_sum = Ten_res2,Bonf_sum=Bonf_sum)
  }
  
  
  return(list(BH_res = BH_res,Bonf_res = Bonf_res,Bonf_sum=Bonf_sum,BH_sum = BH_sum,BH_corr_pvals=BH_corr_pvals))
}

MTR.testversion <- function(p_vals,p_wbh,X_mean2,alpha,p_hat=1,signal=NULL,per_val=0.1){
  
  alpha_hat <- alpha/p_hat
  p_vals2 <- p_vals
  p_vals2[X_mean2<=per_val | X_mean2 >= (1-per_val)] <- 1
  
  W_BH2<- p.adjust(p_wbh,method = "BH")
  T_R <- p_vals2
  R2_BH <- p.adjust(p_vals2[X_mean2>per_val & X_mean2 < (1-per_val)],method = "BH")
  T_R[X_mean2>per_val & X_mean2 < (1-per_val)] <- R2_BH
  R_BH <- p.adjust(p_vals,method = "BH")
  
  BH_corr_pvals <- data.frame(WABH_2 = W_BH2, Ten = T_R,BH = R_BH,ABH = R_BH)
  
  BH_res <- data.frame(WABH_2 = 1*I(W_BH2<alpha_hat), Ten = 1*I(T_R<alpha),BH = 1*I(R_BH<alpha),ABH = 1*I(R_BH<alpha_hat))
  
  W_BY2<- p.adjust(p_wbh,method = "holm")
  R_BY <- p.adjust(p_vals,method = "holm")
  T_R_Bonf <- p_vals2
  R2_BY <- p.adjust(p_vals2[X_mean2>per_val & X_mean2 < (1-per_val)],method = "holm")
  T_R_Bonf[X_mean2>per_val & X_mean2 < (1-per_val)] <- R2_BY
  
  Bonf_res <- data.frame(WABH_2 = 1*I(W_BY2<alpha), Ten = 1*I(T_R_Bonf<alpha),Bonf = 1*I(R_BY<alpha))
  
  BH_sum <- NULL
  Bonf_sum <- NULL
  if(!is.null(signal)){
    M_new <- length(X_mean2)
    n_signal <- !(1:M_new %in% signal)
    BH_sum   <- c(sum(BH_res$BH),sum(BH_res$BH[signal]),sum(BH_res$BH[n_signal]))
    W2_res   <- c(sum(BH_res$WABH_2),sum(BH_res$WABH_2[signal]),sum(BH_res$WABH_2[n_signal]))
    ADAP_res <- c(sum(BH_res$ABH),sum(BH_res$ABH[signal]),sum(BH_res$ABH[n_signal]))
    Ten_res  <- c(sum(BH_res$Ten),sum(BH_res$Ten[signal]),sum(BH_res$Ten[n_signal]))
    
    BH_sum <- data.frame(WABH_sum_2 = W2_res, Ten_sum = Ten_res, ABH_sum = ADAP_res,BH_sum = BH_sum)
    
    W2_res2  <- c(sum(Bonf_res$WABH_2),sum(Bonf_res$WABH_2[signal]),sum(Bonf_res$WABH_2[n_signal]))
    Bonf_sum <- c(sum(Bonf_res$Bonf),sum(Bonf_res$Bonf[signal]),sum(Bonf_res$Bonf[n_signal]))
    Ten_res2 <- c(sum(Bonf_res$Ten),sum(Bonf_res$Ten[signal]),sum(Bonf_res$Ten[n_signal]))
    
    Bonf_sum <- data.frame(WABH_sum_2 = W2_res2, Ten_sum = Ten_res2,Bonf_sum=Bonf_sum)
  }
  
  
  return(list(BH_res = BH_res,Bonf_res = Bonf_res,Bonf_sum=Bonf_sum,BH_sum = BH_sum,BH_corr_pvals=BH_corr_pvals))
}
