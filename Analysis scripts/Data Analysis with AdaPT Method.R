source("Programs/WABH.R")
library(RandomFields)
library(adaptMT)
library(CAMT)
library(ggplot2)
library(qvalue)
library(dplyr)
library(splines)
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


## Removing voxels with damage to less than 0.5% of patients
R_sq <- R_sq[X_mean>=0.005]
lm_info <- lm_info[X_mean>=0.005,]
coord_mat <- coord_mat[X_mean>=0.005,]
X_mean <- X_mean[X_mean>=0.005]

# Study information
M <- dim(lm_info)[1] #number of voxels tested
N <- length(Y_dat) # Sample size
alpha <- 0.05 #Significance level

# Predicted standard error
pred_SE <- 1/sqrt(X_mean*(1-X_mean)*var(Y_dat)*(N-1)/(1-R_sq)) #Sm
pred_SE <- pred_SE*median(lm_info[,2])/median(pred_SE)

#p_values
p_vals <- lm_info[,4]




# Formulas for adapt
X_dat2 <- data.frame(x1 = coord_mat[,1],x2 = coord_mat[,2], x3 = coord_mat[,3], pred_SE = pred_SE)
formulas_mu <- "s(pred_SE)"
formulas_pi <- "s(x1, x2, x3)"
# Running adapt
adapt_test <- adapt_gam(x = X_dat2, pvals = p_vals, pi_formulas = formulas_pi, mu_formulas = formulas_mu, 
                        nfits = 20, alphas = alpha)

# Extracting the estimated non-null probability.
pi_hat <- adapt_test$params[[1]]$pix

# Estimating the WABH weights with the adapt estimate of the non-null probability
tau <- 0.9 # Setting tau, the non-null probability where MMW holds.
wght_padapt_eta <- weighted_p(p_vals,pred_SE, pi_hat, alpha, eta=0.25, p_MMW = tau)
weights_MMW_padapt <- wght_padapt_eta$MMW_weights # Weights

# Getting proportion of inconclusive and upweighted voxels.
length(weights_MMW_padapt[weights_MMW_padapt>1])/length(weights_MMW_padapt)
length(weights_MMW_padapt[weights_MMW_padapt>0.1])/length(weights_MMW_padapt)

# Maximum weight
max(weights_MMW_padapt)

# Performing hypothesis testing on weighted p-values
MTR_testadaptpm_eta <- MTR(p_vals, weights_MMW_padapt$p_weight_MMW, weights_MMW_padapt$p_weight_eta, 
                           X_mean, alpha, mean(1-pi_hat), per_val=0.1)

# Extracting results 
Result_piadapt_MMW <- MTR_testadaptpm_eta$BH_res$WABH

TolRej_piadapt_MMW <- sum(Result_piadapt_MMW)
TolRej_piadapt_MMW

TolRej_adapt <- adapt_test$nrejs
TolRej_adapt

