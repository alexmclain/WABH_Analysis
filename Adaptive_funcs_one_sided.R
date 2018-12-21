

pred_pow_func <- function(g,weights,alpha){
  (1-pnorm(qnorm(1-alpha*weights/length(weights))-g))
}

weight_fun <- function(k,g){
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



### This function computes the FDR constraint in *
FDR.constraint<-function(lambda, g, p = 0.5,alpha = 0.05){
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


### This uses uniroot to find the value of lambda
### that satisfies the constraint above.  May have to mess with
### the endpoints.  Theoretically in (0, infinity).
get.lambda<-function(g, p=.5, alpha=.05){
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

## This puts it all together to return weights

weights_fun_JH<-function(g, p, alpha,find_eta=FALSE){
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
  lambda<-get.lambda(g,p,alpha)
  log_t<- pnorm(log(lambda)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  return(list(log_t=log_t,lambda=lambda))
}

t_eta_finder_jh <- function(s,p,alpha){
  lambda <-get.lambda.eta(s,p,alpha)
  eta <- sqrt(2*log(lambda))*min(s)
  g <- eta/s
  log_t<- pnorm(log(lambda)/g + 0.5*g, lower.tail =FALSE, log.p = TRUE)
  return(list(log_t=log_t,lambda=lambda,eta=eta))
}


### This function computes the FDR constraint in *
FDR.constraint.eta<-function(lambda, s, p = 0.5,alpha = 0.05){
  #This function doesn't require eta, so first eta and "g" are calculated.
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


### This uses uniroot to find the value of lambda
### that satisfies the constraint above.  May have to mess with
### the endpoints.  Theoretically in (0, infinity).
get.lambda.eta<-function(s, p=.5, alpha=.05){
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


#### function to find number of neighbors of TN #####

neigh_func <- function(res,Nei_mat,signal){
  RES <- c(0,0)-1
  M <- length(res)
  res_ind <- 1:M
  rej <- res_ind[res==1]
  TN <- signal[signal %in% rej]
  if(length(TN)>0){
  nei_res <- rep(0,length(TN))
  check <- 1
  for(i in TN){
    nei_res[check] <- sum(1*I(Nei_mat[i,] %in% TN))
    check <- check +1
  }
  RES <- c(mean(nei_res),mean(nei_res>0))
  }
  
  return(RES)
}


Weights <- function(p_vals,eta_pow,pred_SE,p_hat,alpha,plots=FALSE,signal=NULL){
  
  weights <- weight_finder(g = eta_pow/pred_SE,p=p_hat,alpha=alpha)
  p_hat2 <- p_hat
  if(p_hat==1){p_hat2 <- 1-1/M}
  weights2 <- weights_fun_JH(g = eta_pow/pred_SE,p=p_hat2,alpha)
  
  p_weight <- p_vals/weights
  p_weight[p_weight>1] <- 1
  p_weight2 <- p_vals/weights2
  p_weight2[p_weight2>1] <- 1
  
  
  pred_pow3 <- pred_pow_func(g = eta_pow/pred_SE,weights=weights,alpha=alpha)
  pred_pow32<- pred_pow_func(g = eta_pow/pred_SE,weights=weights2,alpha=alpha)
  
  pred_bonf <- sum(pred_pow3*K/M)  #Expected number of Bonferonni rejections. 
  
  if(plots){
    par(mfrow=c(2,2))
    plot(pred_pow3,weights,xlab="Predicted Power",ylab="Weights",ylim=c(0,max(c(weights,weights2))),xlim=range(c(pred_pow3,pred_pow32)))
    points(pred_pow3[signal],weights[signal],col=2,pch=19)
    points(pred_pow32,weights2,col=3,pch=19)
    points(pred_pow32[signal],weights2[signal],col=4,pch=19)
    abline(h=1)
    plot(pred_SE,weights,xlab="Predicted Standard Error",ylab="Weights",ylim=c(0,max(c(weights,weights2))),xlim=c(min(pred_SE),max(pred_SE[is.finite(pred_SE)])))
    points(pred_SE[signal],weights[signal],col=2,pch=19)
    points(pred_SE,weights2,col=3,pch=19)
    points(pred_SE[signal],weights2[signal],col=4,pch=19)
    abline(h=1)
    V1 <- min(c(log(alpha/M),log(p_vals),log(p_weight),log(p_weight2)))
    plot(log(p_vals),ylim=c(V1,0),ylab="Unadjusted P-values",xlab="Position")
    points(signal,log(p_vals[signal]),col=2,pch=19)
    abline(h=log(c(1:5)*alpha/M),lwd=0.5)
    plot(log(p_weight2),ylim=c(V1,0),ylab="Weighted BH P-values",xlab="Position")
    points(signal,log(p_weight2[signal]),col=2,pch=19)
    abline(h=log(c(1:5)*alpha/M),lwd=0.5)
    
  }
  
  return(list(weights_bonf = weights,weights_bh = weights2,p_wbonf=p_weight,p_wbh=p_weight2,pred_bonf=pred_bonf))
}


MTR <- function(p_vals,p_wbonf,p_wbh,X_mean2,alpha,p_hat=1,signal=NULL,per_val=0.1){
  
  alpha_hat <- alpha/p_hat
  p_vals2 <- p_vals
  p_vals2[X_mean2<=per_val] <- 1
  
  W_BH <- p.adjust(p_wbonf,method = "BH")
  W_BH2<- p.adjust(p_wbh,method = "BH")
  T_R <- p_vals2
  R2_BH <- p.adjust(p_vals2[X_mean2>per_val],method = "BH")
  T_R[X_mean2>per_val] <- R2_BH
  R_BH <- p.adjust(p_vals,method = "BH")
  
  BH_corr_pvals <- data.frame(Bonf_BH = W_BH, WABH = W_BH2, Ten = T_R,BH = R_BH,ABH = R_BH)
  
  BH_res <- data.frame(Bonf_BH = 1*I(W_BH<alpha_hat), WABH = 1*I(W_BH2<alpha_hat), Ten = 1*I(T_R<alpha),BH = 1*I(R_BH<alpha),ABH = 1*I(R_BH<alpha_hat))
  
  W_BY <- p.adjust(p_wbonf,method = "holm")
  W_BY2<- p.adjust(p_wbh,method = "holm")
  R_BY <- p.adjust(p_vals,method = "holm")
  T_R_Bonf <- p_vals2
  R2_BY <- p.adjust(p_vals2[X_mean2>per_val],method = "holm")
  T_R_Bonf[X_mean2>per_val] <- R2_BY
  
  Bonf_res <- data.frame(WBonf = 1*I(W_BY<alpha), WABH = 1*I(W_BY2<alpha), Ten = 1*I(T_R_Bonf<alpha),Bonf = 1*I(R_BY<alpha))
  
  BH_sum <- NULL
  Bonf_sum <- NULL
  if(!is.null(signal)){
    n_signal <- !(1:M %in% signal)
    BH_sum   <- c(sum(BH_res$BH),sum(BH_res$BH[signal]),sum(BH_res$BH[n_signal]))
    W_res    <- c(sum(BH_res$Bonf_BH),sum(BH_res$Bonf_BH[signal]),sum(BH_res$Bonf_BH[n_signal]))
    W2_res   <- c(sum(BH_res$WABH),sum(BH_res$WABH[signal]),sum(BH_res$WABH[n_signal]))
    ADAP_res <- c(sum(BH_res$ABH),sum(BH_res$ABH[signal]),sum(BH_res$ABH[n_signal]))
    Ten_res  <- c(sum(BH_res$Ten),sum(BH_res$Ten[signal]),sum(BH_res$Ten[n_signal]))
    
    BH_sum <- data.frame(WBonf_sum = W_res, WABH_sum = W2_res, Ten_sum = Ten_res, ABH_sum = ADAP_res,BH_sum = BH_sum)
    
    W_res2   <- c(sum(Bonf_res$WBonf),sum(Bonf_res$WBonf[signal]),sum(Bonf_res$WBonf[n_signal]))
    W2_res2  <- c(sum(Bonf_res$WABH),sum(Bonf_res$WABH[signal]),sum(Bonf_res$WABH[n_signal]))
    Bonf_sum <- c(sum(Bonf_res$Bonf),sum(Bonf_res$Bonf[signal]),sum(Bonf_res$Bonf[n_signal]))
    Ten_res2 <- c(sum(Bonf_res$Ten),sum(Bonf_res$Ten[signal]),sum(Bonf_res$Ten[n_signal]))
    
    Bonf_sum <- data.frame(WBonf_sum = W_res2, WABH_sum = W2_res2, Ten_sum = Ten_res2,Bonf_sum=Bonf_sum)
  }
  
  
  return(list(BH_res = BH_res,Bonf_res = Bonf_res,Bonf_sum=Bonf_sum,BH_sum = BH_sum,BH_corr_pvals=BH_corr_pvals))
}

