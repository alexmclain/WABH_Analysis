## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

source(paste0("Functions",delim,"WABHProgram.R"))
source(paste0("Functions",delim,"swfdr.R"))
library(adaptMT)
library(RandomFields)
library(ggplot2)
library(cp4p)
library(dplyr)
library(IHW)
library(splines)
library(CAMT)
library(reticulate)

#Reading in different simulation study parameter settings
parlist <- as.matrix(data.frame(read.csv(paste0("Simulation study",delim,"Misspecified Setting Results",delim,"args_list.csv"),header = TRUE, fileEncoding="UTF-8-BOM")))

for(args in 1:54){
  args_list <- as.numeric(unlist(parlist[args,]))
  
  set.seed(66)
  B <- 500 # Number of simulations
  N <- 200 #number of subjects
  x <- seq(1, 100, 1) #Voxel index
  M <- length(x)^2 #Number of voxels
  K <- as.numeric(args_list[1]) #number of false nulls
  #Various data generation parameters#
  err_sd <- 0.8
  alpha <- 0.05
  beta0 <- -1 #alpha0 star in the paper
  eta <- as.numeric(args_list[2])    #theta value in paper draft
  C <- as.numeric(args_list[3]) #Controls the heterogeneity 0.5, 1.5, 3
  sig_sp <- 10 #Spatial clustering of signals
  mis_sp <- as.numeric(args_list[4])*4 #Spatial clustering of misspecified covariate
  lat_sp <- 50 #Spatial clustering of data
  
  
  Weight_res<- matrix(0,B,12)
  Weight_res2<- matrix(0,B,9)
  ihw_mat <- Regular_res <- AD_res <- TenRule_res <- matrix(0,B,3)
  adapt_mat <- swfdr_mat <- camt_mat <- matrix(0,B,3)
  RFoptions(spConform=FALSE)
  start <- as.numeric(args_list[5])
  pi_hat_mean <- NULL
  p3_mean <- NULL
  p_mean <- NULL
  for(k in 1:B){
    set.seed(k + start)
    
    ## Generating the nulls and non-nulls
    signal_data <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp),x=x,y=x,grid=TRUE,n=1)
    quan <- quantile(array(signal_data),prob = 1-K/M)
    signal_data <- 1*I(array(signal_data)>=quan)
    
    ## Generating the latent data (intercepts) alpha0m star in the paper
    data <- RFsimulate(model = RMstable(alpha = 2, scale = lat_sp), x=x, y=x, grid=TRUE,n=1)
    LP_data <- expand.grid(x1=x,x2=x)
    LP_data$Latent_Var <- C*(array(data)-mean(array(data)))/sd(array(data)) #standard deviation is C value
    LP_data$signal <- signal_data+1e-10
    signal<- seq(1,M,1)[signal_data==1]
    
    ## Adding data for misspecification...
    misp_data <- RFsimulate(model = RMstable(alpha = 2, scale = mis_sp),x=x,y=x,grid=TRUE,n=1)
    LP_data$mis_par <- 0.25*(array(misp_data)-mean(array(misp_data)))/sd(array(misp_data))
    mis_cov <- rbinom(N,1,0.5)
    
    ## Generate the X and Y data
    BR_dat <- matrix(0,N,M)
    RE_eff <- err_sd*rnorm(N) # b_i (random effects related to X and Y)
    Y_dat <- 0.5*rnorm(N)+0.5*RE_eff # Y data
    Y_dat <- (Y_dat - mean(Y_dat))/sd(Y_dat)
    eta_samp <- runif(M,0,2*eta) # Alpha coefficients mean is eta
    for(i in 1:M){
      X <- beta0+LP_data$Latent_Var[i]+RE_eff+eta_samp[i]*I(i %in% signal)*Y_dat + LP_data$mis_par[i]*mis_cov
      p_X <- exp(X)/(1+exp(X))
      b_X <- rbinom(N,1,p_X)
      BR_dat[,i] <- b_X
    }
    
    ## X bar m
    X_mean <- apply(BR_dat,2,mean)
    usedtestindicator <- seq(1,M,1)[X_mean*(1-X_mean)>0]
    
    ##Update BR_dat and true value for each test
    BR_dat <- BR_dat[,X_mean*(1-X_mean)>0] 
    signal_data <- signal_data[X_mean*(1-X_mean)>0]
    
    ## Calculating X_plus
    X_pl <- apply(BR_dat,1,sum)
    SIG <- cov(cbind(Y_dat,X_pl))
    
    ## Update X bar m
    X_mean <- apply(BR_dat,2,mean)
    
    ##Creat new M value new test number and update signal value
    
    M_new <- length(X_mean)
    signal<- seq(1,M_new,1)[signal_data==1]
    
    ## Running all the logistic regression models
    R_sq <- rep(0,M_new)
    lm_info <- matrix(0,M_new,4)
    for(i in 1:M_new){
      t_X_pl <- (X_pl-BR_dat[,i])/(M_new-1)
      t_X_pl <- log(t_X_pl/(1-t_X_pl))
      t_X_pl <- (t_X_pl - mean(t_X_pl))/sd(t_X_pl)
      t_lm <- lm(Y_dat~t_X_pl)
      R_sq[i] <- summary(t_lm)$r.squared
      t1 <- try(t_glm <- glm(BR_dat[,i] ~Y_dat+t_X_pl,family=binomial(link='logit')),silent = TRUE)
      if(attr(t1,"class")[1] == "try-error"){t_lm_info <- c(0,1,-999,1)}
      if(attr(t1,"class")[1] == "glm"){t_lm_info <- c(summary(t_glm)$coefficients[2,1:4])}
      lm_info[i,] <- t_lm_info
    }
    
    lm_info[,4] <- pnorm(lm_info[,3],lower.tail = FALSE) #One-sided test
    p_vals <- lm_info[,4]
    p_vals[p_vals==0] <- 1e-20
    p_hat  <- (M_new*pi0est(p_vals, lambda = alpha)$pi0 - 1)/M_new
    
    pred_SE <- 1/sqrt(X_mean*(1-X_mean)*var(Y_dat)*(N-1)/(1-R_sq)) #Sm
    pred_SE <- pred_SE*median(lm_info[,2])/median(pred_SE) 
    X_dat2 <- data.frame(x1 = LP_data$x1[usedtestindicator],x2 = LP_data$x2[usedtestindicator], pred_SE = pred_SE)
    
    ###first set of pim values
    p1 <- rep(1-p_hat,times=M_new)
    p_mean <- c(p_mean, 1-p_hat)
    wght_pconstant <- weighted_p(p_vals,pred_SE, p1, alpha, eta, p_MMW=mean(p1))
    
    # WABH MMW tau constant method Ten Rule method Adaptive BH method and Regular BH method
    MTR_testconstantpm <- MTR(p_vals,wght_pconstant$p_weight_MMW ,wght_pconstant$p_weight_eta,X_mean,alpha,mean(1-p1),signal,per_val=0.1)
    
    
    #Number of discoveries, true discoveries and false discoveries for all procedures
    Weight_res[k,1:3] <- MTR_testconstantpm$BH_sum$WABH_sum # Each row each iteration
    Regular_res[k,] <-  MTR_testconstantpm$BH_sum$BH_sum
    AD_res[k,] <- MTR_testconstantpm$BH_sum$ABH_sum
    TenRule_res[k,] <- MTR_testconstantpm$BH_sum$Ten_sum
    
    ###second set of pim values with AdaPT method
    formulas_mu <- paste0("~ns(pred_SE, df = ", 5, ")")
    formulas_pi <- paste0("~ns(x1, df = ", 5, ")+ns(x2, df = ", 5, ")+ns(x1*x2, df = ", 5, ")")
    
    adapt_test <- adapt_glm(x = X_dat2, pvals = p_vals, pi_formulas = formulas_pi , mu_formulas = formulas_mu, nfits = 5,alphas = alpha)
    
    pi_hat <- adapt_test$params[[1]]$pix
    
    pi_hat_mean <- rbind(pi_hat_mean,c(mean(pi_hat), min(pi_hat), max(pi_hat)) )
    
    wght_padapt <- weighted_p(p_vals,pred_SE, pi_hat, alpha, eta, p_MMW=mean(pi_hat))
    wght_padapt2 <- weighted_p(p_vals,pred_SE, pi_hat, alpha, eta, p_MMW=0.5)
    wght_padapt3 <- weighted_p(p_vals,pred_SE, pi_hat, alpha, eta, p_MMW=0.9)
    
    # WABH MMW tau mean pim method WABH MMW tau 0.5 method
    MTR_testadaptpm <- MTR(p_vals,wght_padapt$p_weight_MMW ,wght_padapt2$p_weight_MMW,X_mean,alpha,mean(1-pi_hat),signal,per_val=0.1)
    
    #Number of discoveries, true discoveries and false discoveries for all procedures
    Weight_res[k,4:6] <-  MTR_testadaptpm$BH_sum$WABH_sum # Each row each iteration
    Weight_res[k,7:9] <-  MTR_testadaptpm$BH_sum$WABH_sum_2 # Each row each iteration
    
    # WABH MMW tau 0.9 method
    MTR_testadaptpm <- MTR(p_vals,wght_padapt3$p_weight_MMW ,wght_padapt2$p_weight_MMW,X_mean,alpha,mean(1-pi_hat),signal,per_val=0.1)
    
    #Number of discoveries, true discoveries and false discoveries for all procedures
    Weight_res[k,10:12] <-  MTR_testadaptpm$BH_sum$WABH_sum # Each row each iteration
    
    ###third set of pim values with CAMT method
    formulas_mu <- paste0("~ns(pred_SE, df = ", 5, ")")
    formulas_pi <- paste0("~ns(x1, df = ", 5, ")+ns(x2, df = ", 5, ")+ns(x1*x2, df = ", 5, ")")
    
    pi0.var <- model.matrix(as.formula(formulas_pi),data=X_dat2)[,-1]
    f1.var <- model.matrix(as.formula(formulas_mu),data=X_dat2)[,-1]
    
    camt_test <- try(camt.fdr(pvals = p_vals, pi0.var = pi0.var, f1.var = f1.var, alg.type = "OS", control.method = "hybrid"), silent = TRUE)
    
    DF <- 5
    while( !is.null(attr(camt_test, "class")) & DF > 0 ){
      DF <- DF - 1
      if(DF > 0){
        formulas_mu <- paste0("~ns(pred_SE, df = ", DF, ")")
        formulas_pi <- paste0("~ns(x1, df = ", DF, ")+ns(x2, df = ", DF, ")+ns(x1*x2, df = ", DF, ")")
        pi0.var <- model.matrix(as.formula(formulas_pi),data=X_dat2)[,-1]
        f1.var <- model.matrix(as.formula(formulas_mu),data=X_dat2)[,-1]
        camt_test <- try(camt.fdr(pvals = p_vals, pi0.var = pi0.var, f1.var = f1.var, alg.type = "OS", control.method = "hybrid"), silent = TRUE)
      }
      
      if(DF == 0){
        pi0.var <- model.matrix(~ pred_SE, data = X_dat2)[,-1]
        f1.var <- model.matrix(~ x1 + x2, data = X_dat2)[,-1]
        camt_test <- try(camt.fdr(pvals = p_vals, pi0.var = pi0.var, f1.var = f1.var, alg.type = "OS", control.method = "hybrid"), silent = TRUE)
      }
      
    }
    
    if(is.null(attr(camt_test, "class")) ){
      p3 <- 1-camt_test$pi0 #third set of pim values with CAMT method
      p3_mean <- rbind(p3_mean, c(mean(p3), min(p3), max(p3)) )
      
      wght_pcamt <- weighted_p(p_vals,pred_SE, p3, alpha, eta, p_MMW=mean(p3))
      wght_pcamt2 <- weighted_p(p_vals,pred_SE, p3, alpha, eta, p_MMW=0.5)
      wght_pcamt3 <- weighted_p(p_vals,pred_SE, p3, alpha, eta, p_MMW=0.9)
      
      
      # WABH MMW tau mean pim method WABH MMW tau 0.5 method
      MTR_testcamt <- MTR(p_vals, wght_pcamt$p_weight_MMW, wght_pcamt2$p_weight_MMW, X_mean, 
                          alpha, mean(1-p3), signal, per_val=0.1)
      
      #Number of discoveries, true discoveries and false discoveries for all procedures
      Weight_res2[k, 1:3] <- MTR_testcamt$BH_sum$WABH_sum # Each row each iteration
      Weight_res2[k, 4:6] <- MTR_testcamt$BH_sum$WABH_sum_2
      
      # WABH MMW tau 0.9 method
      MTR_testcamt <- MTR(p_vals, wght_pcamt3$p_weight_MMW, wght_pcamt2$p_weight_MMW, X_mean, 
                          alpha, mean(1-p3), signal, per_val=0.1)
      
      #Number of discoveries, true discoveries and false discoveries for all procedures
      Weight_res2[k, 7:9] <- MTR_testcamt$BH_sum$WABH_sum
      
      #Number of discoveries, true discoveries and false discoveries for CAMT method
      vals_camt <- I(camt_test$fdr<alpha)*1
      TolRej_camt <- sum(vals_camt)
      if(TolRej_camt>0){
        CorRej_camt <- table(vals_camt,signal_data)[2,2]
        camt_res <- c(TolRej_camt,CorRej_camt,TolRej_camt-CorRej_camt)
        camt_mat[k,] <- camt_res
      }
    }
    
    if(!is.null(attr(camt_test, "class")) ){
      #Setting CAMT results to the adaptive results if it failed.
      Weight_res2[k, 1:3] <- Weight_res[k, 1:3]
      Weight_res2[k, 4:6] <- Weight_res[k, 1:3]
      Weight_res2[k, 7:9] <- Weight_res[k, 1:3]
      camt_mat[k,] <- AD_res[k,]
      cat("CAMT Failed \n")
    }
    
    ###### Fitting the IHW method #####
    ihw_test <- ihw(pvalues=p_vals, covariates=X_dat2$pred_SE, alpha=alpha, covariate_type = "ordinal",nbins = "auto", folds = NULL, quiet = TRUE,nfolds = 5L, nfolds_internal = 5L, nsplits_internal = 1L,lambdas = "auto", seed = 1L, adjustment_type = "BH", null_proportion = TRUE)
    
    TolRej_ihw <- rejections(ihw_test)
    if(TolRej_ihw>0){
      vals_ihw <- rejected_hypotheses(ihw_test)
      CorRej_ihw <- table(vals_ihw,signal_data)[2,2]
      ihw_res <- c(TolRej_ihw,CorRej_ihw,TolRej_ihw-CorRej_ihw)
      ihw_mat[k,] <- ihw_res
    }
    
    ###### Fitting the ADAPT method with GLM #####
    TolRej_adapt <- adapt_test$nrejs
    if(TolRej_adapt>0){
      vals_adapt <- p_vals <= adapt_test$s
      CorRej_adapt <- table(vals_adapt,signal_data)[2,2]
      adapt_res <- c(TolRej_adapt,CorRej_adapt,TolRej_adapt-CorRej_adapt)
      adapt_mat[k,] <- adapt_res
    }
    
    ### Fitting swfdr method ###
    swfdr_test <- lm_qvalue(p_vals, X=X_dat2,type = "linear",smooth.df = 18,smoothing = "smooth.spline")
    
    TolRej_swfdr <- length(p_vals[swfdr_test$qvalues<alpha])
    if(TolRej_swfdr>0){
      vals_swfdr <- I(swfdr_test$qvalues<alpha)*1
      CorRej_swfdr <- table(vals_swfdr,signal_data)[2,2]
      swfdr_res <- c(TolRej_swfdr,CorRej_swfdr,TolRej_swfdr-CorRej_swfdr)
      swfdr_mat[k,] <- swfdr_res
    }
    cat(k,"\n")
  }
  
  res <- rbind( matrix(apply(Weight_res,2,mean),4,3, byrow = TRUE), matrix(apply(Weight_res2,2,mean),3,3, byrow = TRUE), apply(TenRule_res,2,mean), apply(AD_res,2,mean), apply(Regular_res,2,mean),apply(ihw_mat,2,mean),apply(adapt_mat,2,mean),apply(swfdr_mat,2,mean),apply(camt_mat,2,mean))
  
  FDR_vals <- cbind(Weight_res[,3]/Weight_res[,1],Weight_res[,6]/Weight_res[,4],Weight_res[,9]/Weight_res[,7],Weight_res[,12]/Weight_res[,10],Weight_res2[,3]/Weight_res2[,1],Weight_res2[,6]/Weight_res2[,4],Weight_res2[,9]/Weight_res2[,7],TenRule_res[,3]/TenRule_res[,1],AD_res[,3]/AD_res[,1],Regular_res[,3]/Regular_res[,1],ihw_mat[,3]/ihw_mat[,1],adapt_mat[,3]/adapt_mat[,1],swfdr_mat[,3]/swfdr_mat[,1],camt_mat[,3]/camt_mat[,1])
  
  #The FDR value should be 0 when the total rejection is 0
  FDR_vals[Weight_res[,1]==0,1] <- 0
  FDR_vals[Weight_res[,4]==0,2] <- 0
  FDR_vals[Weight_res[,7]==0,3] <- 0
  FDR_vals[Weight_res[,10]==0,4] <- 0
  FDR_vals[Weight_res2[,1]==0,5] <- 0
  FDR_vals[Weight_res2[,4]==0,6] <- 0
  FDR_vals[Weight_res2[,7]==0,7] <- 0
  FDR_vals[TenRule_res[,1]==0,8] <- 0
  FDR_vals[AD_res[,1]==0,9] <- 0
  FDR_vals[Regular_res[,1]==0,10] <- 0
  FDR_vals[ihw_mat[,1]==0,11] <- 0
  FDR_vals[adapt_mat[,1]==0,12] <- 0
  FDR_vals[swfdr_mat[,1]==0,13] <- 0
  FDR_vals[camt_mat[,1]==0,14] <- 0

  
  res <- cbind(res,apply(FDR_vals,2,mean))
  
  rownames(res) <- c("WABH MMW constant pm","WABH MMW mean_pm adapt pm","WABH MMW 0.5 adapt pm","WABH MMW 0.9 adapt pm","WABH MMW mean_pm camt pm","WABH MMW 0.5 camt pm","WABH MMW 0.9 camt pm","Ten Rule","Adaptive BH","Regular BH","IHW","ADAPT","SWFDR","CAMT")
  colnames(res) <- c("Total Rej","Correct Rej","False Rej","FDR")
  res
  
  res_full <- cbind(Weight_res,Weight_res2,TenRule_res,AD_res,Regular_res,ihw_mat,adapt_mat,swfdr_mat,camt_mat)
  
  
  c_names <- c("WABH MMW constant pm","WABH MMW mean_pm adapt pm","WABH MMW 0.5 adapt pm","WABH MMW 0.9 adapt pm","WABH MMW mean_pm camt pm","WABH MMW 0.5 camt pm","WABH MMW 0.9 camt pm","Ten Rule","Adaptive BH","Regular BH","IHW","ADAPT","SWFDR","CAMT")
  colnames(res_full) <- paste0(rep(c_names,each = 3), rep(c("_Tot","_Cor","_False"),length(c_names)))
  rownames(res_full) <- NULL
  res_full
  
  write.csv(res,file=paste0("Simulation study",delim,"Misspecified Setting Results",delim,"Summarized",delim,"Res_mis K=",K," eta=",eta," C=",C," mis_sp=",mis_sp," start=",start,".csv"))
  write.csv(res_full,file=paste0("Simulation study",delim,"Misspecified Setting Results",delim,"By_iteration",delim,"Res_full_mis K=",K," eta=",eta," C=",C," mis_sp=",mis_sp," start=",start,".csv"))
}

