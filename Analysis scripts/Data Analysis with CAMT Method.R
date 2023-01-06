source("WABH.R")
library(RandomFields)
library(ggplot2)
library(qvalue)
library(dplyr)
library(splines)
library(CAMT)
library(mgcv)

# Outcome data
Y_dat <- readRDS("Data/Y_dat.rds")

# Reading in the coordinates and results from the logistic model.
coord_mat_reg_results <- readRDS("Data/Log_mod_res_logitY_coord_mat.rds")
head(coord_mat_reg_results)

# Extracting key data elements.
coord_mat <- coord_mat_reg_results[,8:10] #3D Coordinates
X_mean <- coord_mat_reg_results[,2] # Proportion of damaged lesions
R_sq <- coord_mat_reg_results[,3]  # R-squared
lm_info <- coord_mat_reg_results[,4:7] #Logistic regression information (coef est, SE, Z, p-value)


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


# Formulas for CAMT
X_dat2 <- data.frame(x1 = coord_mat[,1],x2 = coord_mat[,2], x3 = coord_mat[,3], pred_SE = pred_SE)
DF <- 12
formulas_mu <- paste0("~ns(pred_SE, df = ", DF, ")")
formulas_pi <- paste0("~ns(x1, df = ", DF, ")+ns(x2, df = ", DF, ")+ns(x3, df = ", DF, ")+ns(x1*x2, df = ", DF, ")+ns(x1*x3, df = ", DF, ")+ns(x2*x3, df = ", DF, ")")

pi0.var <- model.matrix(as.formula(formulas_pi),data=X_dat2)[,-1]
f1.var <- model.matrix(as.formula(formulas_mu),data=X_dat2)[,-1]

camt_test <- camt.fdr(pvals = p_vals, pi0.var = pi0.var, f1.var = f1.var, 
                      alg.type = "EM", control.method = "hybrid", pi0.low = 0.1)

# Extracting the estimated non-null probability.
p3 <- 1-camt_test$pi0


tau <- 0.9 # Setting tau, the non-null probability where MMW holds.
wght_camt_eta <- weighted_p(p_vals,pred_SE, p3, alpha, eta=0.25, p_MMW = tau)
weights_MMW_pcamt <- wght_camt_eta$MMW_weights # Weights

# Getting proportion of inconclusive and upweighted voxels for MMW.
length(weights_MMW_pcamt[weights_MMW_pcamt>1])/length(weights_MMW_pcamt)
length(weights_MMW_pcamt[weights_MMW_pcamt>0.1])/length(weights_MMW_pcamt)

# Maximum weight
max(weights_MMW_pcamt)

# Performing hypothesis testing on weighted p-values
MTR_testcamtpm_eta <- MTR(p_vals, wght_camt_eta$p_weight_MMW, wght_camt_eta$p_weight_eta, 
                          X_mean, alpha, mean(1-p3), per_val=0.1)

#Number of discoveries for all procedures
Result_picamt_MMW <- MTR_testcamtpm_eta$BH_res$WABH

TolRej_picamt_MMW <- sum(Result_picamt_MMW)
TolRej_picamt_MMW

vals_camt <- I(camt_test$fdr<alpha)*1
TolRej_camt <- sum(vals_camt)
TolRej_camt



#10% Rule weights calculation
ind_10per <- I(X_mean>=0.1 & X_mean <=0.9)*1
weights_10per <- (M*ind_10per)/sum(ind_10per)


# Getting proportion of inconclusive and upweighted voxels for 10% rule.
length(weights_10per[weights_10per>1])/length(weights_10per)
length(weights_10per[weights_10per>0.1])/length(weights_10per)

# Max weight for 10% rule.
max(weights_10per)

Result_Ten <- MTR_testcamtpm_eta0.25$BH_res$Ten
Result_BH <- MTR_testcamtpm_eta0.25$BH_res$BH
Result_ABH <- MTR_testcamtpm_eta0.25$BH_res$ABH

TolRej_Ten <- sum(Result_Ten)
TolRej_BH <- sum(Result_BH)
TolRej_ABH <- sum(Result_ABH)





