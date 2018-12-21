
source("Adaptive_funcs_one_sided.R")
library(dplyr)
library(qvalue)

link_data <- read.csv("id_vec.csv")
link_data <- link_data[,-1]
link_data2 <- link_data[link_data$NoDup,]

Cov_dat <- read.csv("ID_w_Y_Xp_nodup.csv")
Cov_dat <- Cov_dat[,-1]
#cov1 = Months post stroke at date of examination.
#cov2 = Age at test.
Cov_dat <- Cov_dat[,-5] #Remove Age at test to limit missing data.

#Use the id-vec to properly order the data.
CD = merge(link_data2,Cov_dat,by.x = "x",by.y = "ID",all.x = TRUE,sort = FALSE)
CD20 <- CD[complete.cases(CD),]

alpha <- 0.1
M <- 936457
Y_dat <- CD20$Y_dat
cov1 <- CD20$cov1
X_pl <- CD20$Xp
ID_NA <- CD20$x

Y_dat <- (100-Y_dat)/100
Y_dat <- log(Y_dat/(1-Y_dat))
N <- length(Y_dat)

con = file("Small_long_img.txt", "r")

#Run the analysis and save the results (might take a while).

X_mean <- rep(0,M)
R_sq <- rep(0,M)
lm_info <- matrix(0,M,4)
for(i in 1:M){
  t1 <- c(as.double(read.table(con, nrows=1)))
  BR_dat <- t1[link_data$NoDup]
  BR_dat <- BR_dat[complete.cases(CD)]
  X_mean[i] <- mean(BR_dat)
  if(X_mean[i]==0){
    t_lm_info <- c(0,1,0,1)
    cat(i,"zero mean error \n")}
  if(X_mean[i]>0){
    t_X_pl <- (X_pl-BR_dat)/(M-1)
    t_X_pl <- log(t_X_pl/(1-t_X_pl))
    t_X_pl <- (t_X_pl - mean(t_X_pl))/sd(t_X_pl)
    t1 <- try(t_lm <- lm(Y_dat~t_X_pl+cov1))
    if(attr(t1,"class")[1] == "try-error"){
      t_lm_info <- c(0,1,0,1)
      cat(i,"LM error \n")}
    if(attr(t1,"class")[1] == "lm"){
      R_sq[i] <- summary(t_lm)$r.squared
      t1 <- try(t_glm <- glm(BR_dat ~Y_dat+t_X_pl+cov1,family=binomial(link='logit')),silent = TRUE)
      if(attr(t1,"class")[1] == "try-error"){
        t_lm_info <- c(0,1,0,1)
        cat(i,"GLM error \n")}
      if(attr(t1,"class")[1] == "glm"){
        t_lm_info <- c(summary(t_glm)$coefficients[2,1:4])}
    }
  }
  lm_info[i,] <- t_lm_info
}



lm_info[,4] <- pnorm(lm_info[,3],lower.tail = FALSE)

write.csv(cbind(X_mean,R_sq,lm_info),file = "Log_mod_res_logitY_w_MSS.csv")

close(con)

dat <- read.csv("Log_mod_res_logitY_w_MSS.csv")
X_mean <- dat[,2]
R_sq <- dat[,3]
lm_info <- dat[,4:7]


P_lm_info <- lm_info

X_mean <- X_mean[lm_info[,1]!=0]
R_sq <- R_sq[lm_info[,1]!=0]
lm_info <- lm_info[lm_info[,1]!=0,]

M <- dim(lm_info)[1]

### Calculate the predicted SE ###
pred_SE <- 1/sqrt(X_mean*(1-X_mean)*var(Y_dat)*(N-1)/(1-R_sq))
pred_SE[X_mean==0] <- 1

#### Estimate the null proportion of p-values ####
p_vals <- lm_info[,4]
p_hat  <- pi0est(p_vals,lambda = alpha)$pi0
p      <- p_hat

##### Finding the zero derivative value of eta #####
weights1 <- weights_fun_JH(g = pred_SE,p=p,alpha=alpha,find_eta = TRUE)

eta_pow<- weights1$eta
weights2<- weights1$weights

#eta_pow <- 0.25
#weights <- weights_fun_JH(g = eta_pow/pred_SE,p=p,alpha=alpha)$weights

## Showing the the weights are ordered WRT their SE.
s_w = weights2[order(pred_SE)]
diff = s_w[-length(s_w)] - s_w[-1]
#plot(s_w,type="l")
any(diff<0)

## Get the weighted p_values
p_wbh <- p_weight2 <- p_vals/weights2
p_wbh[p_wbh>1] <- 1
p_weight2[p_weight2>1] <- 1

## Get the multiple testing results
MTR_res <- MTR(p_vals,p_wbonf,p_wbh,X_mean,alpha,p_hat,per_val = 0.1)

## Summarize the results 
apply(MTR_res$BH_res,2,sum)
apply(MTR_res$Bonf_res,2,sum)

## Get the predicted number of Bonferonni rejections
## Estimated number of true nulls
K <- M*p_hat


#pred_pow3 <- pred_pow_func(g = eta_pow/pred_SE,weights=weights,alpha=alpha)
# Get the power  for each test
pred_pow32<- pred_pow_func(g = eta_pow/pred_SE,weights=weights2,alpha=alpha)
#Expected number of Bonferonni rejections. 
sum(pred_pow32*(1-K/M))
#Proportion of weights greater then 1.
length(weights2[weights2>1])/length(weights2)

# Put the results together
fin_res <- cbind(P_lm_info,0,0,0,0,0,1,1,1)
fin_res[P_lm_info[,1]!=0,4] <- p_vals
fin_res[P_lm_info[,1]!=0,5] <- weights2
fin_res[P_lm_info[,1]!=0,6:10] <- MTR_res$BH_res
fin_res[P_lm_info[,1]!=0,11] <- p_weight2
fin_res[P_lm_info[,1]!=0,12] <- pvals_10

#Summary of weights for significant voxels
summary(weights2)
summary(weights2[MTR_res$BH_corr_pvals$WABH*p_hat<alpha])

fwrite(as.list(fin_res),"Final_initial_res.csv")
  




