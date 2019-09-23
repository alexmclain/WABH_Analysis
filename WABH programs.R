


weights_fun_JH<-function(g, p, alpha,find_eta=FALSE){
  ### This function implements the optimal weight estimation discussed in Section 3.2 of the main paper.
  ### Arguments:
  # g = eta_m/S_m if find_eta=FALSE, or
  #   = S_m       if find_eta=TRUE.  Here, S_m is the predicted standard error and eta_m the beta value.
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  # find_eta: set to true if the MMW criteria should be applied and eta be calculated 
  #           (in this case g=S_m), set to false if the eta is given and g = eta_m/S_m.
  ### Value:
  # The function will return the estimated weights, the value of lambda (c in equation (4)) and the   # value of eta (this is noteworthy if find_eta=TRUE).
  
  
  if(!find_eta){
    t_find <-t_finder_jh(g,p,alpha)
    thresholds <- exp(t_find$log_t)
    eta=NULL
    lambda = t_find$lambda
  }
  if(find_eta){
    t_eta_find <-t_eta_finder_jh(g,p,alpha)
    thresholds <- exp(t_eta_find$log_t)
    eta=t_eta_find$eta   
    lambda = t_eta_find$lambda 
  }
  return(list(weights = thresholds/mean(thresholds),lambda=lambda,eta=eta))
}


t_finder_jh <- function(g,p,alpha){
  ## This function calculate the value of lambda (c in equation (4)) and then calculates the
  ## weights.  For this function eta is given.
  ### Arguments:
  # g = eta_m/S_m.  Here, S_m is the predicted standard error and eta_m the beta value.
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  
  lambda<-get.lambda(g,p,alpha)  ### Getting the lambda for the Lagrange multiplier
  log_t<- pnorm(log(lambda)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE) ## Weights
  return(list(log_t=log_t,lambda=lambda))
}

t_eta_finder_jh <- function(s,p,alpha){
  ## This function calculate the value of eta that satisfies the MMW criteria, lambda (c in
  ## equation (4)) and then calculates the weights.  
  ### Arguments:
  # g = S_m.  Here, S_m is the predicted standard error.
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  
  lambda <-get.lambda.eta(s,p,alpha) ### Getting the lambda for the Lagrange multiplier
  eta <- sqrt(2*log(lambda))*min(s)  ### Value of eta that satisfies MMW.
  g <- eta/s
  log_t<- pnorm(log(lambda)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE) ## Weights
  return(list(log_t=log_t,lambda=lambda,eta=eta))
}



get.lambda<-function(g, p=.5, alpha=.05){
  ### This uses uniroot to find the value of lambda
  ### that satisfies the constraint above.  May have to mess with
  ### the endpoints.  Theoretically in (0, infinity).
  ### Arguments:
  # g: eta/S_m.  
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  FDR.c<-function(lambda)
  {
    FDR.constraint(lambda,g,p, alpha)
  }
  t1 <- try(asd <- uniroot(FDR.c, c(1e-100, length(g)/alpha)),silent = TRUE)
  if(class(t1)!="try-error"){lambda <<- asd$root}
  if(class(t1)=="try-error"){
    t1 <- try(asd <- uniroot(FDR.c, c(1e-100, 500)),silent = TRUE)
    if(class(t1)!="try-error"){lambda <<- asd$root}
    if(class(t1)=="try-error"){stop(t1[1])}
  }
  return(lambda)
}


FDR.constraint<-function(lambda, g, p = 0.5,alpha = 0.05){
  ### This function computes the FDR constraint in given in equation (4) for t defined in equation (3).
  ### Arguments:
  # lambda: Lagrange multiplier (c in equation (4) of the paper.)
  # g: eta/S_m.  
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  
  #Calculate the threshold on the log scale
  log_t<- pnorm(log(lambda)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  log_t[log_t>0] <- 0
  t <- exp(log_t)
  #Calculate the log-power then the power.
  log_pow <- pnorm(qnorm(log_t,lower.tail = FALSE, log.p = TRUE)-g, lower.tail =FALSE, log.p = TRUE)
  pow <- exp(log_pow)
  #Calculate the expected FDR.
  efdr <- p*mean(t)/(p*mean(t) + (1-p)*mean(pow))
  return(log(efdr) - log(alpha))
}


get.lambda.eta<-function(s, p=.5, alpha=.05){
  ### This uses uniroot to find the value of lambda for the MMW criteria
  ### that satisfies the FDR constraint above.  May have to mess with
  ### the endpoints.  Theoretically in (0, infinity).
  ### Arguments:
  # s = S_m.  Here, S_m is the predicted standard error.
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  FDR.c<-function(lambda)
  {
    FDR.constraint.eta(lambda,s,p, alpha)
  }
  t1 <- try(asd <- uniroot(FDR.c, c(1.01, length(s)/alpha)),silent = TRUE)
  if(class(t1)!="try-error"){lambda <<- asd$root}
  if(class(t1)=="try-error"){
    t1 <- try(asd <- uniroot(FDR.c, c(1e-100, 500)),silent = TRUE)
    if(class(t1)!="try-error"){lambda <<- asd$root}
    if(class(t1)=="try-error"){stop(t1[1])}
  }
  return(lambda)
}



pred_pow_func <- function(g,weights,alpha){
  # predicted power for an effect size of g and alpha level = alpha*weights/length(weights).
  (1-pnorm(qnorm(1-alpha*weights/length(weights))-g))
}

weight_fun <- function(k,g){
  # equation (3) in the paper.
  # g = eta_m/S_m
  # k = log(c)
  w_m <- pnorm(g*0.5 + k/g,lower.tail = FALSE,log.p = TRUE)
  w_m[w_m>0] <- 0
  exp(w_m)
}

weight_max <- function(k,g,p,alpha){
  
  (alpha/p - sum(weight_fun(k,g)))
}

weight_finder <- function(g,p,alpha=0.05){
  t <- t_finder_acm(g,p,alpha)
  return(t/mean(t))
}

t_finder_acm <- function(g,p,alpha){
  asd <- uniroot(weight_max, c(1e-100, 99999999999999),g=g,p=p,alpha=alpha)
  k_ast <- asd$root
  t <- weight_fun(k_ast,g)
  t
}




## This puts it all together to return weights


FDR.constraint.eta<-function(lambda, s, p = 0.5,alpha = 0.05){
  ### This function computes the FDR constraint in given in equation (4) for t defined in equation (5).
  ### Arguments:
  # lambda: Lagrange multiplier (c in equation (4) of the paper.)
  # s: S_m.  
  # p: estimated proportion of null hypotheses.
  # alpha: significance level.
  
  #This function doesn't require eta, so first eta and "g" are calculated for the given lambda.
  eta  <- sqrt(2*log(lambda))*min(s)
  g    <- eta/s
  
  #Calculate the threshold on the log scale
  log_t<- pnorm(log(lambda)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  log_t[log_t>0] <- 0
  t <- exp(log_t)
  #Calculate the log-power then the power.
  log_pow <- pnorm(qnorm(log_t,lower.tail = FALSE, log.p = TRUE)-g, lower.tail =FALSE, log.p = TRUE)
  pow <- exp(log_pow)
  #Calculate the expected FDR.
  efdr <- p*mean(t)/(p*mean(t) + (1-p)*mean(pow))
  return(log(efdr) - log(alpha))
}




MTR <- function(p_vals,p_wbonf,p_wbh,X_mean2,alpha,p_hat=1,signal=NULL,per_val=0.1){
  ### This function calculates the results of many multiple testing procedures WABH, BH, ABH and the p% rule
  ### Arguements:
  # p_vals: unadjusted p-values
  # p_wbonf: weighted p-values
  # p_wbh: weighted p-values
  # X_mean2: mean of the data (used for the p% rule)
  # alpha: significance level
  # p_hat: estimated proportion of null hypotheses
  # signal: which tests are actually signals (used in simulations)
  # per_val: at what proportion the the p% rule cutoff (default is per_val=0.1 for the 10% rule)
  ### Values:
  # This function returns a list which contains:
  #       - BH_res: full results from the BH procedure where 0 indicates non-significant and 1 
  #         indicates significant.  The columns of the results are for the WABH with p_wbonf, 
  #         WABH with p_wbh, the p% rule, the standard unweighted BH and the standard unweighted ABH.
  #       - Bonf_res: full results from the Holm procedure where 0 indicates non-significant and 1 
  #         indicates significant.  The columns of the results are for p_wbonf, p_wbh, the p% rule
  #         and unweighted (no adative procedures).
  #       - BH_sum and Bonf_sum (only if signal is given) give the number of true and false discoveries.
  #       - BH_corr_pvals: the BH adjusted p-values for p_wbonf, p_wbh, the p% rule, and standard BH.
  #         Note that these adjusted p-values to not take the estimated proportion of null into account.
  
  
  alpha_hat <- alpha/p_hat
  p_vals2 <- p_vals
  p_vals2[X_mean2<=per_val] <- 1
  
  W_BH <- p.adjust(p_wbonf,method = "BH")
  W_BH2<- p.adjust(p_wbh,method = "BH")
  T_R <- p_vals2
  R2_BH <- p.adjust(p_vals2[X_mean2>per_val],method = "BH")
  T_R[X_mean2>per_val] <- R2_BH
  R_BH <- p.adjust(p_vals,method = "BH")
  
  BH_corr_pvals <- data.frame(WABH = W_BH, WABH_2 = W_BH2, Ten = T_R,BH = R_BH,ABH = R_BH)
  
  BH_res <- data.frame(WABH = 1*I(W_BH<alpha_hat), WABH_2 = 1*I(W_BH2<alpha_hat), Ten = 1*I(T_R<alpha),BH = 1*I(R_BH<alpha),ABH = 1*I(R_BH<alpha_hat))
  
  W_BY <- p.adjust(p_wbonf,method = "holm")
  W_BY2<- p.adjust(p_wbh,method = "holm")
  R_BY <- p.adjust(p_vals,method = "holm")
  T_R_Bonf <- p_vals2
  R2_BY <- p.adjust(p_vals2[X_mean2>per_val],method = "holm")
  T_R_Bonf[X_mean2>per_val] <- R2_BY
  
  Bonf_res <- data.frame(WABH = 1*I(W_BY<alpha), WABH_2 = 1*I(W_BY2<alpha), Ten = 1*I(T_R_Bonf<alpha),Bonf = 1*I(R_BY<alpha))
  
  BH_sum <- NULL
  Bonf_sum <- NULL
  if(!is.null(signal)){
    n_signal <- !(1:M %in% signal)
    BH_sum   <- c(sum(BH_res$BH),sum(BH_res$BH[signal]),sum(BH_res$BH[n_signal]))
    W_res    <- c(sum(BH_res$Bonf_BH),sum(BH_res$Bonf_BH[signal]),sum(BH_res$Bonf_BH[n_signal]))
    W2_res   <- c(sum(BH_res$WABH),sum(BH_res$WABH[signal]),sum(BH_res$WABH[n_signal]))
    ADAP_res <- c(sum(BH_res$ABH),sum(BH_res$ABH[signal]),sum(BH_res$ABH[n_signal]))
    Ten_res  <- c(sum(BH_res$Ten),sum(BH_res$Ten[signal]),sum(BH_res$Ten[n_signal]))
    
    BH_sum <- data.frame(WABH_sum = W_res, WABH_sum_2 = W2_res, Ten_sum = Ten_res, ABH_sum = ADAP_res,BH_sum = BH_sum)
    
    W_res2   <- c(sum(Bonf_res$WBonf),sum(Bonf_res$WBonf[signal]),sum(Bonf_res$WBonf[n_signal]))
    W2_res2  <- c(sum(Bonf_res$WABH),sum(Bonf_res$WABH[signal]),sum(Bonf_res$WABH[n_signal]))
    Bonf_sum <- c(sum(Bonf_res$Bonf),sum(Bonf_res$Bonf[signal]),sum(Bonf_res$Bonf[n_signal]))
    Ten_res2 <- c(sum(Bonf_res$Ten),sum(Bonf_res$Ten[signal]),sum(Bonf_res$Ten[n_signal]))
    
    Bonf_sum <- data.frame(WABH_sum = W_res2, WABH_sum_2 = W2_res2, Ten_sum = Ten_res2,Bonf_sum=Bonf_sum)
  }
  return(list(BH_res = BH_res,Bonf_res = Bonf_res,Bonf_sum=Bonf_sum,BH_sum = BH_sum,BH_corr_pvals=BH_corr_pvals))
}

