## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

library(RandomFields)
library(ggplot2)
library(cp4p)
library(dplyr)
library(splines)
library(reticulate)

parameterlist <- data.frame(seed = 1:5, C = c(3,0.5,1.5,rep(1.5,2)), s = c(rep(10,3),0.01,5) )


L <- nrow(parameterlist)
for (j in 1:L){
args_list <- as.numeric(unlist(parameterlist[j,]))

set.seed(66)
B <- 500 # Number of simulations
N <- 200 #number of subjects
x <- seq(1, 100, 1) #Voxel index
M <- length(x)^2 #Number of voxels
K <- 500 #number of false nulls
#Various data generation parameters#
err_sd <- 0.8
alpha <- 0.05
beta0 <- -1 #alpha0 star in the paper
eta <- 0.25    #theta value in paper draft
C <- as.numeric(args_list[2]) #Controls the heterogeneity 0.5, 1.5, 3
sig_sp <- as.numeric(args_list[3]) #Spatial clustering of signals
lat_sp <- 50 #Spatial clustering of data

simulationno <- as.numeric(args_list[1])
RFoptions(spConform=FALSE)

set.seed(simulationno)

## Generating the nulls and non-nulls
signal_data <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp),x=x,y=x,grid=TRUE,n=1)
quan <- quantile(array(signal_data),prob = 1-K/M)
signal_data <- 1*I(array(signal_data)>=quan)

## Generating the latent data (intercepts) alpha0m star in the paper
data <- RFsimulate(model = RMstable(alpha = 2, scale = lat_sp), x=x, y=x, grid=TRUE,n=1)
LP_data <- expand.grid(x1=x,x2=x)
LP_data$Latent_Var <- C*(array(data)-mean(array(data)))/sd(array(data)) + beta0 #standard deviation is C value
LP_data$signal <- signal_data 
signal<- seq(1,M,1)[signal_data==1]

## Generate the X and Y data
BR_dat <- matrix(0,N,M)
RE_eff <- err_sd*rnorm(N) # b_i (random effects related to X and Y)
Y_dat <- 0.5*rnorm(N)+0.5*RE_eff # Y data
Y_dat <- (Y_dat - mean(Y_dat))/sd(Y_dat)
eta_samp <- runif(M,0,2*eta) # Alpha coefficients mean is eta
for(i in 1:M){
  X <- LP_data$Latent_Var[i]+RE_eff+eta_samp*I(i %in% signal)*Y_dat
  p_X <- exp(X)/(1+exp(X))
  b_X <- rbinom(N,1,p_X)
  BR_dat[,i] <- b_X
}

LP_data$lesionstatus <- BR_dat[100,]
LP_data$signal1[LP_data$signal==0] <- 0
LP_data$signal1[LP_data$signal==1] <- 1

pdf(paste0("Simulation study",delim,"Figures",delim,"Data_examples",delim,"Non Null ",simulationno,"K500theta",eta,"C",C,"sig_sp",sig_sp,".pdf"),height = 6.75,width = 5)


print(ggplot(data = LP_data, aes(x = x1, y = x2, fill = signal1)) + geom_tile() + 
  
  scale_fill_gradientn(colours = c("gray90","black"), values = c(0,1)) + 
  
  xlab(expression(z[m1]^pi)) + ylab(expression(z[m2]^pi)) +
  
  theme(axis.title.x = element_text(size=18, face="bold"),
        
        axis.title.y = element_text(size=18, face="bold"),axis.text.x = element_text(size=14),
        
        axis.text.y = element_text(size=14), panel.background = element_rect(fill = "white", colour = "grey50"),
        
        legend.position="none"))

dev.off()

pdf(paste0("Simulation study",delim,"Figures",delim,"Data_examples",delim,"Intercept ",simulationno,"K500theta",eta,"C",C,"sig_sp",sig_sp,".pdf"),height = 6.75,width = 6)

print(ggplot(data = LP_data, aes(x = x1, y = x2)) + geom_tile(aes(fill = Latent_Var)) +
  
        xlab(expression(z[m1]^pi)) + ylab(expression(z[m2]^pi)) +
  
  theme(axis.title.x = element_text(size=18, face="bold"),
        
        axis.title.y = element_text(size=18, face="bold"),axis.text.x = element_text(size=14),
        
        axis.text.y = element_text(size=14), panel.background = element_rect(fill = "white", colour = "grey50"),
        
        legend.position="right") + guides(fill=guide_legend(title=expression(alpha[0]^"*"))))

dev.off()


pdf(paste0("Simulation study",delim,"Figures",delim,"Data_examples",delim,"X value ",simulationno,"K500theta",eta,"C",C,"sig_sp",sig_sp,".pdf"),height = 6.75,width = 5)

print(ggplot(data = LP_data, aes(x = x1, y = x2)) + geom_tile(aes(fill = lesionstatus)) +
        
        scale_fill_gradientn(colours = c("gray90","black"), values = c(0,1)) + 
        
        xlab(expression(z[m1]^pi)) + ylab(expression(z[m2]^pi))  +
        
        theme(axis.title.x = element_text(size=18, face="bold"),
              
              axis.title.y = element_text(size=18, face="bold"),axis.text.x = element_text(size=14),
              
              axis.text.y = element_text(size=14), panel.background = element_rect(fill = "white", colour = "grey50"),
              
              legend.position="none"))

dev.off()

}
