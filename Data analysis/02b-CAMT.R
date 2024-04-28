## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

source(paste0("Functions",delim,"WABHProgram.R"))
library(CAMT)
library(ggplot2)
library(qvalue)
library(dplyr)
library(splines)
library(mgcv)
library(data.table)

# Outcome data
Y_dat <- readRDS(paste0("Data analysis",delim,"Data",delim,"Y_dat.rds"))

# Reading in the coordinates and results from the logistic model.
coord_mat_reg_results <- readRDS(paste0("Data analysis",delim,"Results",delim,"Log_mod_res_logitY_coord_mat.rds"))
head(coord_mat_reg_results)

# Extracting key data elements.
coord_mat <- coord_mat_reg_results[,7:9] #3D Coordinates
X_mean <- coord_mat_reg_results[,1] # Proportion of damaged lesions
R_sq <- coord_mat_reg_results[,2]  # R-squared
lm_info <- coord_mat_reg_results[,3:6] #Logistic regression information (coef est, SE, Z, p-value)


# Removing all voxels with no damage
R_sq <- R_sq[X_mean!=0]
lm_info <- lm_info[X_mean!=0,]
coord_mat <- coord_mat[X_mean!=0,]
X_mean <- X_mean[X_mean!=0]

# Study information
M <- dim(lm_info)[1] #number of voxels tested
N <- length(Y_dat) # Sample size
alpha <- 0.05 #Significance level

# Predicted standard error
pred_SE <- 1/sqrt(X_mean*(1-X_mean)*var(Y_dat)*(N-1)/(1-R_sq)) #Sm
pred_SE <- pred_SE*median(lm_info[,2])/median(pred_SE)

#p_values
p_vals <- lm_info[,4]

remove( list = c("lm_info","coord_mat_reg_results") )
gc()


# Formulas for CAMT
X_dat2 <- data.frame(x1 = coord_mat[,1],x2 = coord_mat[,2], x3 = coord_mat[,3], pred_SE = pred_SE)
DF <- 12
formulas_mu <- paste0("~ns(pred_SE, df = ", DF, ")")
formulas_pi <- paste0("~ns(x1, df = ", DF, ")+ns(x2, df = ", DF, ")+ns(x3, df = ", DF, ")+ns(x1*x2, df = ", DF, ")+ns(x1*x3, df = ", DF, ")+ns(x2*x3, df = ", DF, ")")

pi0.var <- model.matrix(as.formula(formulas_pi),data=X_dat2)[,-1]
f1.var <- model.matrix(as.formula(formulas_mu),data=X_dat2)[,-1]

# The camt.fdr function will take ~ 12 hours to run.
# Below will upload the results of this function as a .rds file.
camt_test <- camt.fdr(pvals = p_vals, pi0.var = pi0.var, f1.var = f1.var, 
                      alg.type = "EM", control.method = "hybrid", pi0.low = 0.1,
                      EM.paras = list(iterlim = 200))

camt_test <- readRDS(paste0("Data analysis",delim,"Results",delim,"CAMT_results.rds"))

# Extracting the estimated non-null probability.
pi_camt <- 1-camt_test$pi0

# WABH CAMT MMW tau 0.9 weights and number of discoveries
# Estimating the WABH weights with the CAMT estimate of the non-null probability
tau <- max(pi_camt) # Setting tau, the non-null probability where MMW holds.
wght_camt_eta <- weighted_p(p_vals,pred_SE, pi_camt, alpha, eta=0.25, p_MMW = tau)
weights_MMW_pcamt <- wght_camt_eta$MMW_weights # WABH CAMT MMW tau 0.9 Weights

# Getting proportion of inconclusive and upweighted voxels.
wgtpro4camttau0.9 <- length(weights_MMW_pcamt[weights_MMW_pcamt>4])/length(weights_MMW_pcamt)
wgtpro2camttau0.9 <- length(weights_MMW_pcamt[weights_MMW_pcamt>2])/length(weights_MMW_pcamt)
wgtpro1camttau0.9 <- length(weights_MMW_pcamt[weights_MMW_pcamt>1])/length(weights_MMW_pcamt)
wgtpro0.1camttau0.9 <- length(weights_MMW_pcamt[weights_MMW_pcamt>0.1])/length(weights_MMW_pcamt)

# Maximum weight
wgtmaxcamttau0.9 <- max(weights_MMW_pcamt)

# Performing hypothesis testing on weighted p-values
MTR_testcamtpm_etatau0.9 <- MTR(p_vals, wght_camt_eta$p_weight_MMW, wght_camt_eta$p_weight_eta, 
                          X_mean, alpha, mean(1-pi_camt), per_val=0.1)

# Extracting results

Result_picamt_MMW <- MTR_testcamtpm_etatau0.9$BH_res$WABH

TolRej_picamt_MMWtau0.9 <- sum(Result_picamt_MMW)

# WABH CAMT MMW tau 0.5 weights and number of discoveries
# Estimating the WABH weights with the CAMT estimate of the non-null probability
tau1 <- 0.5 # Setting tau, the non-null probability where MMW holds.
wght_camt_etatau0.5 <- weighted_p(p_vals,pred_SE, pi_camt, alpha, eta=0.25, p_MMW = tau1)
weights_MMWtau0.5_pcamt <- wght_camt_etatau0.5$MMW_weights # WABH CAMT MMW tau 0.5 Weights

# Getting proportion of inconclusive and upweighted voxels.
wgtpro4camttau0.5 <- length(weights_MMWtau0.5_pcamt[weights_MMWtau0.5_pcamt>4])/length(weights_MMWtau0.5_pcamt)
wgtpro2camttau0.5 <- length(weights_MMWtau0.5_pcamt[weights_MMWtau0.5_pcamt>2])/length(weights_MMWtau0.5_pcamt)
wgtpro1camttau0.5 <- length(weights_MMWtau0.5_pcamt[weights_MMWtau0.5_pcamt>1])/length(weights_MMWtau0.5_pcamt)
wgtpro0.1camttau0.5 <- length(weights_MMWtau0.5_pcamt[weights_MMWtau0.5_pcamt>0.1])/length(weights_MMWtau0.5_pcamt)

# Maximum weight
wgtmaxcamttau0.5 <- max(weights_MMWtau0.5_pcamt)

# Performing hypothesis testing on weighted p-values
MTR_testcamtpm_etatau0.5 <- MTR(p_vals, wght_camt_etatau0.5$p_weight_MMW, wght_camt_etatau0.5$p_weight_eta, 
                                X_mean, alpha, mean(1-pi_camt), per_val=0.1)

# Extracting results

Result_picamt_MMWtau0.5 <- MTR_testcamtpm_etatau0.5$BH_res$WABH

TolRej_picamt_MMWtau0.5 <- sum(Result_picamt_MMWtau0.5)

finalresultscamt <- matrix(c(wgtpro4camttau0.5,wgtpro2camttau0.5,wgtpro1camttau0.5,wgtpro0.1camttau0.5,wgtmaxcamttau0.5,TolRej_picamt_MMWtau0.5, 
                             wgtpro4camttau0.9,wgtpro2camttau0.9,wgtpro1camttau0.9,wgtpro0.1camttau0.9,wgtmaxcamttau0.9,TolRej_picamt_MMWtau0.9),ncol=6,byrow=TRUE)

finalresultscamt <- as.data.frame(round(finalresultscamt,3))
rownames(finalresultscamt) <- c('WABH MMW tau0.5 CAMT','WABH MMW tau0.9 CAMT')
colnames(finalresultscamt) <- c('proportion of weights >= 4','proportion of weights >= 2','proportion of weights >= 1','proportion of weights >= 0.1','maximum weight','number of discoveries')


## Results in Table 1 of supplemental material (last 2 rows).
finalresultscamt

## Percent of discoveries with weights greater than 4
table(weights_MMW_pcamt[Result_picamt_MMW == 1]>4)/length(weights_MMW_pcamt[Result_picamt_MMW == 1]>4)

#Number of discoveries for CAMT
vals_camt <- I(camt_test$fdr<alpha)*1
TolRej_camt <- sum(vals_camt)
TolRej_camt

#10% Rule weights calculation
ind_10per <- I(X_mean>=0.1 & X_mean <=0.9)*1
weights_10per <- (M*ind_10per)/sum(ind_10per)


#Number of discoveries for 10% rule, BH and ABH method
(TolRej_Ten <- sum(MTR_testcamtpm_etatau0.9$BH_res$Ten))
(TolRej_BH <- sum(MTR_testcamtpm_etatau0.9$BH_res$BH))
(TolRej_ABH <- sum(MTR_testcamtpm_etatau0.9$BH_res$ABH))



### Exporting results ###
coord_mat_reg_results <- readRDS(paste0("Data analysis",delim,"Results",delim,"Log_mod_res_logitY_coord_mat.rds"))
lm_info <- coord_mat_reg_results[,3:6] #Logistic regression information (coef est, SE, Z, p-value)
X_mean <- coord_mat_reg_results[,1] # Proportion of damaged lesions
lm_info <- lm_info[X_mean!=0,]
X_mean <- X_mean[X_mean!=0]

fin_res <- cbind(lm_info,0,0,0,0,0,0,0,1,0,0,0)
fin_res[,4] <- p_vals
fin_res[,5] <- weights_MMW_pcamt
fin_res[,6] <- Result_picamt_MMW
fin_res[,7] <- Result_picamt_MMWtau0.5
fin_res[,8] <- Result_Ten
fin_res[,9] <- Result_BH
fin_res[,10] <- Result_ABH
fin_res[,11] <- vals_camt
fin_res[,12] <- wght_camt_eta$p_weight_MMW
fin_res[,13] <- pi_camt
fin_res[,14] <- camt_test$k
fin_res[,15] <- weights_MMWtau0.5_pcamt


saveRDS(fin_res,paste0("Data analysis",delim,"Results",delim,"Final_initial_res.rds"))


