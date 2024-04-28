## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

source(paste0("Functions",delim,"WABHProgram.R"))
library(adaptMT)
library(ggplot2)
library(qvalue)
library(dplyr)
library(splines)
library(mgcv)

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

remove( list = c("lm_info","coord_mat_reg_results") )
gc()

# Formulas for adapt
X_dat2 <- data.frame(x1 = coord_mat[,1],x2 = coord_mat[,2], x3 = coord_mat[,3], pred_SE = pred_SE)
formulas_mu <- "s(pred_SE)"
formulas_pi <- "s(x1, x2, x3)"
# Running adapt

# The adapt_gam function will take ~ 3 weeks to run.
# Below will upload the results of this function as a .rds file.
adapt_test <- adapt_gam(x = X_dat2, pvals = p_vals, pi_formulas = formulas_pi, mu_formulas = formulas_mu, 
                        nfits = 20, alphas = alpha, verbose = list(print = TRUE, fit = TRUE, ms = TRUE))

adapt_test <- readRDS(paste0("Data analysis",delim,"Results",delim,"Adapt_results.rds"))

# AdaPT number of discoveries
TolRej_adapt <- adapt_test$nrejs
TolRej_adapt

# Extracting the estimated non-null probability.
pi_hat <- adapt_test$params[[1]]$pix

# WABH AdaPT MMW tau 0.9 weights and number of discoveries
# Estimating the WABH weights with the adapt estimate of the non-null probability
tau <- max(pi_hat) # Setting tau, the non-null probability where MMW holds.
wght_padapt_eta <- weighted_p(p_vals,pred_SE, pi_hat, alpha, eta=0.25, p_MMW = tau)
weights_MMW_padapt <- wght_padapt_eta$MMW_weights # WABH AdaPT MMW tau 0.9 Weights

# Getting proportion of inconclusive and upweighted voxels.
wgtpro4adapttau0.9 <- length(weights_MMW_padapt[weights_MMW_padapt>4])/length(weights_MMW_padapt)
wgtpro2adapttau0.9 <- length(weights_MMW_padapt[weights_MMW_padapt>2])/length(weights_MMW_padapt)
wgtpro1adapttau0.9 <- length(weights_MMW_padapt[weights_MMW_padapt>1])/length(weights_MMW_padapt)
wgtpro0.1adapttau0.9 <- length(weights_MMW_padapt[weights_MMW_padapt>0.1])/length(weights_MMW_padapt)

# Maximum weight
wgtmaxadapttau0.9 <- max(weights_MMW_padapt)

# Performing hypothesis testing on weighted p-values
MTR_testadaptpm_eta <- MTR(p_vals, wght_padapt_eta$p_weight_MMW, wght_padapt_eta$p_weight_eta, 
                           X_mean, alpha, mean(1-pi_hat), per_val=0.1)

# Extracting results 
Result_piadapt_MMW <- MTR_testadaptpm_eta$BH_res$WABH

TolRej_piadapt_MMWtau0.9 <- sum(Result_piadapt_MMW)


# WABH AdaPT MMW tau 0.5 weights and number of discoveries
# Estimating the WABH weights with the adapt estimate of the non-null probability
tau1 <- 0.5 # Setting tau, the non-null probability where MMW holds.
wght_padapt_etatau0.5 <- weighted_p(p_vals,pred_SE, pi_hat, alpha, eta=0.25, p_MMW = tau1)
weights_MMWtau0.5_padapt <- wght_padapt_etatau0.5$MMW_weights # WABH AdaPT MMW tau 0.5 Weights

# Getting proportion of inconclusive and upweighted voxels.
wgtpro4adapttau0.5 <- length(weights_MMWtau0.5_padapt[weights_MMWtau0.5_padapt>4])/length(weights_MMWtau0.5_padapt)
wgtpro2adapttau0.5 <- length(weights_MMWtau0.5_padapt[weights_MMWtau0.5_padapt>2])/length(weights_MMWtau0.5_padapt)
wgtpro1adapttau0.5 <- length(weights_MMWtau0.5_padapt[weights_MMWtau0.5_padapt>1])/length(weights_MMWtau0.5_padapt)
wgtpro0.1adapttau0.5 <- length(weights_MMWtau0.5_padapt[weights_MMWtau0.5_padapt>0.1])/length(weights_MMWtau0.5_padapt)

# Maximum weight
wgtmaxadapttau0.5 <- max(weights_MMWtau0.5_padapt)

# Performing hypothesis testing on weighted p-values
MTR_testadaptpm_etatau0.5 <- MTR(p_vals, wght_padapt_etatau0.5$p_weight_MMW, wght_padapt_etatau0.5$p_weight_eta, 
                           X_mean, alpha, mean(1-pi_hat), per_val=0.1)

# Extracting results 
Result_piadapt_MMWtau0.5 <- MTR_testadaptpm_etatau0.5$BH_res$WABH

TolRej_piadapt_MMWtau0.5 <- sum(Result_piadapt_MMWtau0.5)


#Truncating pi_hat at 0.9
pi_hat_0.9cut <- pi_hat
pi_hat_0.9cut[pi_hat_0.9cut>0.9] <- 0.9

# WABH AdaPT0.9cut MMW tau 0.9 weights and number of discoveries
# Estimating the WABH weights with the adapt estimate of the non-null probability
tau <- max(pi_hat_0.9cut) # Setting tau, the non-null probability where MMW holds.
wght_padapt0.9cut_etatau0.9 <- weighted_p(p_vals,pred_SE, pi_hat_0.9cut, alpha, eta=0.25, p_MMW = tau)
weights_MMWtau0.9_padapt0.9cut <- wght_padapt0.9cut_etatau0.9$MMW_weights # WABH AdaPT0.9cut MMW tau 0.9 Weights

# Getting proportion of inconclusive and upweighted voxels.
wgtpro4adapt0.9cuttau0.9 <- length(weights_MMWtau0.9_padapt0.9cut[weights_MMWtau0.9_padapt0.9cut>4])/length(weights_MMWtau0.9_padapt0.9cut)
wgtpro2adapt0.9cuttau0.9 <- length(weights_MMWtau0.9_padapt0.9cut[weights_MMWtau0.9_padapt0.9cut>2])/length(weights_MMWtau0.9_padapt0.9cut)
wgtpro1adapt0.9cuttau0.9 <- length(weights_MMWtau0.9_padapt0.9cut[weights_MMWtau0.9_padapt0.9cut>1])/length(weights_MMWtau0.9_padapt0.9cut)
wgtpro0.1adapt0.9cuttau0.9 <- length(weights_MMWtau0.9_padapt0.9cut[weights_MMWtau0.9_padapt0.9cut>0.1])/length(weights_MMWtau0.9_padapt0.9cut)

# Maximum weight
wgtmaxadapt0.9cuttau0.9 <- max(weights_MMWtau0.9_padapt0.9cut)

# Performing hypothesis testing on weighted p-values
MTR_testadapt0.9cutpm_etatau0.9 <- MTR(p_vals, wght_padapt0.9cut_etatau0.9$p_weight_MMW, wght_padapt0.9cut_etatau0.9$p_weight_eta, 
                           X_mean, alpha, mean(1-pi_hat_0.9cut), per_val=0.1)

# Extracting results 
Result_piadapt0.9cut_MMWtau0.9 <- MTR_testadapt0.9cutpm_etatau0.9$BH_res$WABH

TolRej_piadapt0.9cut_MMWtau0.9 <- sum(Result_piadapt0.9cut_MMWtau0.9)

# WABH AdaPT0.9cut MMW tau 0.5 weights and number of discoveries
# Estimating the WABH weights with the adapt estimate of the non-null probability
tau1 <- 0.5 # Setting tau, the non-null probability where MMW holds.
wght_padapt0.9cut_etatau0.5 <- weighted_p(p_vals,pred_SE, pi_hat_0.9cut, alpha, eta=0.25, p_MMW = tau1)
weights_MMWtau0.5_padapt0.9cut <- wght_padapt0.9cut_etatau0.5$MMW_weights # WABH AdaPT0.9cut MMW tau 0.5 Weights

# Getting proportion of inconclusive and upweighted voxels.
wgtpro4adapt0.9cuttau0.5 <- length(weights_MMWtau0.5_padapt0.9cut[weights_MMWtau0.5_padapt0.9cut>4])/length(weights_MMWtau0.5_padapt0.9cut)
wgtpro2adapt0.9cuttau0.5 <- length(weights_MMWtau0.5_padapt0.9cut[weights_MMWtau0.5_padapt0.9cut>2])/length(weights_MMWtau0.5_padapt0.9cut)
wgtpro1adapt0.9cuttau0.5 <- length(weights_MMWtau0.5_padapt0.9cut[weights_MMWtau0.5_padapt0.9cut>1])/length(weights_MMWtau0.5_padapt0.9cut)
wgtpro0.1adapt0.9cuttau0.5 <- length(weights_MMWtau0.5_padapt0.9cut[weights_MMWtau0.5_padapt0.9cut>0.1])/length(weights_MMWtau0.5_padapt0.9cut)

# Maximum weight
wgtmaxadapt0.9cuttau0.5 <- max(weights_MMWtau0.5_padapt0.9cut)

# Performing hypothesis testing on weighted p-values
MTR_testadapt0.9cutpm_etatau0.5 <- MTR(p_vals, wght_padapt0.9cut_etatau0.5$p_weight_MMW, wght_padapt0.9cut_etatau0.5$p_weight_eta, 
                                       X_mean, alpha, mean(1-pi_hat_0.9cut), per_val=0.1)

# Extracting results 
Result_piadapt0.9cut_MMWtau0.5 <- MTR_testadapt0.9cutpm_etatau0.5$BH_res$WABH

TolRej_piadapt0.9cut_MMWtau0.5 <- sum(Result_piadapt0.9cut_MMWtau0.5)

finalresults <- matrix(c(wgtpro4adapttau0.5,wgtpro2adapttau0.5,wgtpro1adapttau0.5,wgtpro0.1adapttau0.5,wgtmaxadapttau0.5,TolRej_piadapt_MMWtau0.5,
                         wgtpro4adapttau0.9,wgtpro2adapttau0.9,wgtpro1adapttau0.9,wgtpro0.1adapttau0.9,wgtmaxadapttau0.9,TolRej_piadapt_MMWtau0.9,
                         wgtpro4adapt0.9cuttau0.5,wgtpro2adapt0.9cuttau0.5,wgtpro1adapt0.9cuttau0.5,wgtpro0.1adapt0.9cuttau0.5,wgtmaxadapt0.9cuttau0.5,TolRej_piadapt0.9cut_MMWtau0.5,
                         wgtpro4adapt0.9cuttau0.9,wgtpro2adapt0.9cuttau0.9,wgtpro1adapt0.9cuttau0.9,wgtpro0.1adapt0.9cuttau0.9,wgtmaxadapt0.9cuttau0.9,TolRej_piadapt0.9cut_MMWtau0.9),ncol=6,byrow=TRUE)

finalresults <- as.data.frame(round(finalresults,3))
rownames(finalresults) <- c('WABH MMW tau0.5 AdaPT','WABH MMW tau0.9 AdaPT','WABH MMW tau0.5 AdaPT0.9cut','WABH MMW tau0.9 AdaPT0.9cut')
colnames(finalresults) <- c('proportion of weights >= 4','proportion of weights >= 2','proportion of weights >= 1','proportion of weights >= 0.1','maximum weight','number of discoveries')

## Results in Table 1 of supplemental material (first 4 rows).
finalresults

## Percent of discoveries with weights greater than 4
table(weights_MMWtau0.9_padapt0.9cut[Result_piadapt0.9cut_MMWtau0.9 == 1]>4)/length(weights_MMWtau0.9_padapt0.9cut[Result_piadapt0.9cut_MMWtau0.9 == 1]>4)


# Number of discoveries for AdaPT
TolRej_adapt <- adapt_test$nrejs
TolRej_adapt
TolRej_adapt/length(adapt_test$qvals)
