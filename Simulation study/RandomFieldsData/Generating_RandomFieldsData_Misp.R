## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

setwd("C:/Users/mclaina/OneDrive - University of South Carolina/Research/Imaging/Multiple testing/Programs/Adaptive Power/WABH_Analysis")

delim <- "/"

source(paste0("Functions",delim,"WABHProgram.R"))
source(paste0("Functions",delim,"swfdr.R"))
source(paste0("Functions",delim,"law_funcs.R"))
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
  RFoptions(spConform=FALSE)
  start <- as.numeric(args_list[5])
  
  misp_data_mat <- NULL
  
  for(k in 1:B){
    set.seed(k + start)
    
    ## Generating the nulls and non-nulls
    signal_data <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp),x=x,y=x,grid=TRUE,n=1)
    quan <- quantile(array(signal_data),prob = 1-K/M)
    signal_data <- 1*I(array(signal_data)>=quan)
    signal<- seq(1,M,1)[signal_data==1]
    
    ## Generating the latent data (intercepts) alpha0m star in the paper
    data <- RFsimulate(model = RMstable(alpha = 2, scale = lat_sp), x=x, y=x, grid=TRUE,n=1)
    LP_data <- expand.grid(x1=x,x2=x)
    LP_data$Latent_Var <- beta0+C*(array(data)-mean(array(data)))/sd(array(data)) #standard deviation is C value
    LP_data$signal <- signal_data+1e-10
    
    ## Generating data for misspecified settings...
    misp_data <- RFsimulate(model = RMstable(alpha = 2, scale = mis_sp),x=x,y=x,grid=TRUE,n=1)
    misp_data_mat <- rbind(misp_data_mat,  round(array(misp_data), 2) )
  }
  
  saveRDS(misp_data_mat,paste0("Simulation study",delim,"RandomFieldsData",delim,"misp_data K=",K," eta=",eta," C=",C," mis_sp=",mis_sp," start=",start,".rds"))
  cat(args, "\n")
}

