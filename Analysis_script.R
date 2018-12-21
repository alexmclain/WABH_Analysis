
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

#What value of eta to choose (take a look at the applied paper)?
#Plot the significant voxels for Naive, 10% rule and proposed method.
#Plot the power, for the choosen value of eta. This plot is crucial for all analyses that use our method.

P_lm_info <- lm_info

X_mean <- X_mean[lm_info[,1]!=0]
R_sq <- R_sq[lm_info[,1]!=0]
lm_info <- lm_info[lm_info[,1]!=0,]

M <- dim(lm_info)[1]


pred_SE <- 1/sqrt(X_mean*(1-X_mean)*var(Y_dat)*(N-1)/(1-R_sq))
pred_SE[X_mean==0] <- 1

length(X_mean[X_mean>0.1])

ter <- lm(pred_SE~lm_info[,2])
summary(ter)

p_vals <- lm_info[,4]
p_hat  <- pi0est(p_vals,lambda = alpha)$pi0
p      <- p_hat

##### Finding the zero derivative value of eta #####
weights1 <- weights_fun_JH(g = pred_SE,p=p,alpha=alpha,find_eta = TRUE)

eta_pow<- weights1$eta
weights2<- weights1$weights

weights <- weight_finder(g = eta_pow/pred_SE,p=p,alpha=alpha)

eta_pow <- 0.25
weights4 <- weights_fun_JH(g = eta_pow/pred_SE,p=p,alpha=alpha)$weights

s_w = weights2[order(pred_SE)]
diff = s_w[-length(s_w)] - s_w[-1]
#plot(s_w,type="l")
any(diff<0)


p_wbonf <- p_weight <- p_vals/weights
p_wbh <- p_weight2 <- p_vals/weights2
p_wbonf[p_wbonf>1] <- 1
p_wbh[p_wbh>1] <- 1
p_weight[p_weight>1] <- 1
p_weight2[p_weight2>1] <- 1

MTR_res <- MTR(p_vals,p_wbonf,p_wbh,X_mean,alpha,p_hat,per_val = 0.2)
apply(MTR_res$BH_res,2,sum)
apply(MTR_res$Bonf_res,2,sum)


K <- M*p_hat


pred_pow3 <- pred_pow_func(g = eta_pow/pred_SE,weights=weights,alpha=alpha)
pred_pow32<- pred_pow_func(g = eta_pow/pred_SE,weights=weights2,alpha=alpha)

sum(pred_pow3*(1-K/M))  #Expected number of Bonferonni rejections. 
sum(pred_pow32*(1-K/M))
length(weights2[weights2>1])/length(weights2)

weights4 <- weights2

  pdf("Pow_plot_logit.pdf")
  
  par(mfrow=c(1,2))
  plot(pred_pow3[order(pred_pow3,decreasing = FALSE)],weights[order(pred_pow3,decreasing = FALSE)],type="l",xlab="Predicted Power",ylab="Weights",ylim=c(0,max(c(weights,weights2))),xlim=range(c(pred_pow3,pred_pow32)))
  lines(pred_pow32[order(pred_pow32,decreasing = FALSE)],weights2[order(pred_pow32,decreasing = FALSE)],col=3)
  abline(h=1)
  
  plot(pred_SE[order(pred_SE,decreasing = FALSE)],weights[order(pred_SE,decreasing = FALSE)],type = "l",xlab="Predicted Standard Error",ylab="Weights",ylim=c(0,max(c(weights,weights2))),xlim=c(min(pred_SE),max(pred_SE[is.finite(pred_SE)])))
  lines(pred_SE[order(pred_SE,decreasing = FALSE)],weights2[order(pred_SE,decreasing = FALSE)],col=3)
  abline(h=1)
  V1 <- min(c(log(alpha/M),log(p_vals),log(p_weight),log(p_weight2)))
  
  dev.off()
  
  Weight_10 <- rep(length(weights)/length(X_mean[X_mean>0.1]),length(weights))
  Weight_10[X_mean<=0.1] <- 0
  
  Weight_20 <- rep(length(weights)/length(X_mean[X_mean>0.2]),length(weights))
  Weight_20[X_mean<=0.2] <- 0
  
  pred_pow33<- pred_pow_func(g = eta_pow/pred_SE,weights=Weight_10,alpha=alpha)
  sum(pred_pow33*(1-K/M))
  length(Weight_10[Weight_10>1])/length(Weight_10)
  pred_pow33<- pred_pow_func(g = eta_pow/pred_SE,weights=Weight_20,alpha=alpha)
  sum(pred_pow33*(1-K/M))
  length(Weight_20[Weight_20>1])/length(Weight_20)
  
  
  pdf("Weight_plot_logit_w_MPS.pdf",width = 11,height = 6)
  
  par(mfrow=c(1,1))
  plot(pred_SE[order(pred_SE,decreasing = FALSE)],weights2[order(pred_SE,decreasing = FALSE)],type = "l",xlab="Predicted Standard Error",ylab="Weights",lwd=3,xlim=c(min(pred_SE),0.15),ylim=c(0,max(weights3)),cex.lab=1.4,cex.axis = 1.3)
  lines(pred_SE[order(pred_SE,decreasing = FALSE)],Weight_10[order(pred_SE,decreasing = FALSE)],lty=1,col="gray70",lwd=3)
#  lines(pred_SE[order(pred_SE,decreasing = FALSE)],Weight_20[order(pred_SE,decreasing = FALSE)],lty=1,col="gray70",lwd=3)
  lines(pred_SE[order(pred_SE,decreasing = FALSE)],weights3[order(pred_SE,decreasing = FALSE)],lty=3,lwd=3)
  lines(pred_SE[order(pred_SE,decreasing = FALSE)],weights4[order(pred_SE,decreasing = FALSE)],lty=2,lwd=3)
  abline(h=1)
  #legend(0.12,10,legend = c(expression(eta == 0.3), expression(eta == 0.17),expression(eta == 0.1),"10% Rule","20% Rule"),lwd=3,col=c(2,1,3:5),cex=1.5)
  
  dev.off()
  
  
  
  V1 <- min(c(log(alpha/M,base=10),log(p_vals,base=10),log(p_weight2,base=10)))
  pdf("Manhattan_plot_logit_w_MPS.pdf",width = 12,height = 8)
  
  par(mfrow=c(1,2))
  plot(-log(p_vals,base=10),ylim=c(0,-V1),xlim=c(1,M),ylab=expression(-log[10](P)),xlab="Position",pch='.',cex=0.1)
  abline(h=-log(c(1)*alpha/M,base=10),lwd=2)
  
  plot(-log(p_weight2,base=10),ylim=c(0,-V1),ylab=expression(-log[10](Q)),xlab="Position",pch='.',cex=0.1)
  abline(h=-log(c(1)*alpha/M,base=10),lwd=2)
  
  dev.off()
  
  pvals_10 <- p_vals/Weight_10
  pvals_20 <- p_vals/Weight_20
  pvals_10[pvals_10>1] <- 1
  pvals_20[pvals_20>1] <- 1
  
  pdf("Manhattan_plot_logit_w_Age_1020rules.pdf",width = 8,height = 8)
  
  par(mfrow=c(1,2))
  plot(-log(pvals_10),ylim=c(0,-V1),xlim=c(1,M),ylab=expression(-log(p)),xlab="Position",pch='.',cex=0.1)
  abline(h=-log(c(1)*alpha/M),lwd=2)
  
  plot(-log(pvals_20),ylim=c(0,-V1),ylab=expression(-log(p)),xlab="Position",pch='.',cex=0.1)
  abline(h=-log(c(1)*alpha/M),lwd=2)
  
  dev.off()
  
  

    
fin_res <- cbind(P_lm_info,0,0,0,0,0,1,1,1)
fin_res[P_lm_info[,1]!=0,4] <- p_vals
fin_res[P_lm_info[,1]!=0,5] <- weights2
fin_res[P_lm_info[,1]!=0,6:10] <- MTR_res$BH_res
fin_res[P_lm_info[,1]!=0,11] <- p_weight2
fin_res[P_lm_info[,1]!=0,12] <- pvals_10

summary(weights2[MTR_res$BH_corr_pvals$WABH*p_hat<alpha])
summary(weights2[MTR_res$BH_corr_pvals$WABH*p_hat<alpha]<=1)

summary(weights2[MTR_res$BH_corr_pvals$WABH*p_hat<alpha]<=min(weights2[weights2<Weight_20]))
summary(weights2[MTR_res$BH_corr_pvals$WABH*p_hat<alpha]>=max(Weight_20))
summary(weights2>=max(Weight_20))

fwrite(as.list(fin_res),"Final_initial_res.csv")
  




