source("WABH programs.R")
library(dplyr)
library(qvalue)

#Loading in voxel-level mean (proportion of lesions), R-sq, and logistic model information (estimated regression coefficients, standard error, z-values and p-values).

reg_results <- read.csv("Log_mod_res_logitY_w_MSS.csv")

X_mean <- reg_results[,2] # Proportion of lesions
R_sq <- reg_results[,3]  # R-squared
lm_info <- reg_results[,4:7] #Logistic regression information

P_lm_info <- lm_info


# Removing all voxels with no damage
R_sq <- R_sq[X_mean!=0]
lm_info <- lm_info[X_mean!=0,]
X_mean <- X_mean[X_mean!=0]

M <- dim(lm_info)[1] #number of voxels tested
N <- 182 # Sample size
alpha <- 0.1 #Significance level


#predicting standard error, using equation (1) in the paper, where 3.597663 is the variance of the Y data.
pred_SE <- 1/sqrt(X_mean*(1-X_mean)*3.597663*(N-1)/(1-R_sq))

#Getting the number 
#length(X_mean[X_mean>0.1])

#Getting the p_values and estimating the proportion of nulls.
p_vals <- lm_info[,4]
p_hat  <- pi0est(p_vals,lambda = alpha)$pi0
p      <- p_hat


##### Estiamting the p-value weights using the methods discussed in Section 3.2 #####
##### Estimating the weights using the monotone minimum weight (MMW) criteria in Theorem 1 #####
MMW_weights <- weights_fun_JH(g = pred_SE,p=p,alpha=alpha,find_eta = TRUE)

eta_pow_MMW<- MMW_weights$eta #the value of eta that satisfys the MMW criteria 
weights_MMW <- MMW_weights$weights #MMW weights



##### Estimating the weights using specified values of eta: #####
eta_pow <- 0.25
weights_0.25 <- weights_fun_JH(g = eta_pow/pred_SE,p=p,alpha=alpha)$weights
eta_pow <- 0.1
weights_0.1 <- weights_fun_JH(g = eta_pow/pred_SE,p=p,alpha=alpha)$weights

##### The weights for the "10% Rule" ######
Weight_10 <- rep(length(weights_MMW)/length(X_mean[X_mean>0.1]),length(weights_MMW))
Weight_10[X_mean<=0.1] <- 0

### Creating Figure 1 in the paper ###
par(mfrow=c(1,1),mar = c(5,5,2,2))
plot(pred_SE[order(pred_SE,decreasing = FALSE)],weights_MMW[order(pred_SE,decreasing = FALSE)],type = "l",xlab=expression(paste("Predicted Standard Error ",(S[m]))),ylab=expression(paste("Weights ",(w[m]))),lwd=3,xlim=c(min(pred_SE),0.15),ylim=c(0,max(weights_0.1)),cex.lab=1.4,cex.axis = 1.3)
lines(pred_SE[order(pred_SE,decreasing = FALSE)],Weight_10[order(pred_SE,decreasing = FALSE)],lty=1,col="gray70",lwd=3)
lines(pred_SE[order(pred_SE,decreasing = FALSE)],weights_0.1[order(pred_SE,decreasing = FALSE)],lty=3,lwd=3)
lines(pred_SE[order(pred_SE,decreasing = FALSE)],weights_0.25[order(pred_SE,decreasing = FALSE)],lty=2,lwd=3)
abline(h=1)



##### Estimating the proportion of positive weights, max weight and expected number of Bonferroni rejections (all given in table 1)

length(weights_MMW[weights_MMW>1])/length(weights_MMW)
length(weights_0.25[weights_0.25>1])/length(weights_0.25)
length(Weight_10[Weight_10>1])/length(Weight_10)

max(weights_MMW)
max(weights_0.25)
max(Weight_10)

pred_pow_MMW <- pred_pow_func(g = eta_pow_MMW/pred_SE,weights=weights_MMW,alpha=alpha)
pred_pow_0.25 <- pred_pow_func(g = 0.25/pred_SE,weights=weights_0.25,alpha=alpha)
K <- M*p_hat

sum(pred_pow_MMW*(1-K/M))  #Expected number of Bonferonni rejections. 
sum(pred_pow_0.25*(1-K/M))



## Now to calculating the results of the multiple testing.  The function MTR will calculate the results for the BH and adative BH procedure for two weighted p-values and for the 10% rule (using per_val = 0.1). It also calculates the results for the Holm procedure which we'll not use here.

# Calculating the weighted p-values
p_wbh_MMW <- p_vals/weights_MMW

p_wbh_0.25 <- p_vals/weights_0.25

MTR_res <- MTR(p_vals,p_wbh_MMW,p_wbh_0.25,X_mean,alpha,p_hat,per_val = 0.1)
BH_summary <- apply(MTR_res$BH_res,2,sum)
## This gives the number of rejections for the MMW, eta=0.25, the 10% rule, along with the standard (unweighted) BH and adaptive BH procedures
BH_summary

## Getting the some summary statistics about the weights for the significant voxels.
summary(weights_MMW[MTR_res$BH_corr_pvals$WABH*p_hat<alpha])
summary(weights_MMW[MTR_res$BH_corr_pvals$WABH*p_hat<alpha]>3) # Number > 3
mean(weights_MMW[MTR_res$BH_corr_pvals$WABH*p_hat<alpha]>3)  # Proportion > 3


#### Merging the results from voxels with non-zero proportion of damage to the entire 
#### dataset so the results can be exporting and ploted on a brain map.

fin_res <- cbind(P_lm_info,0,0,0,0,0,1,1,1)
fin_res[P_lm_info[,1]!=0,4] <- p_vals
fin_res[P_lm_info[,1]!=0,5] <- weights2
fin_res[P_lm_info[,1]!=0,6:10] <- MTR_res$BH_res
fin_res[P_lm_info[,1]!=0,11] <- p_wbh_MMW
fin_res[P_lm_info[,1]!=0,12] <- p_wbh_0.25

fwrite(as.list(fin_res),"Final_initial_res.csv")


