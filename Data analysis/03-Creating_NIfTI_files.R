## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

library(neurobase)
library(oro.nifti)

fin_res <- readRDS(paste0("Data analysis",delim,"Results",delim,"Final_initial_res.rds"))
big_to_small_data <- readRDS(paste0("Data analysis",delim,"Data",delim,"big_to_small_data.rds"))
table(big_to_small_data)


tot_res <- as.matrix(cbind(big_to_small_data,0,0,0,0,0,0,0,0,0,0,0,0,0))
new_dat <- as.matrix(cbind(fin_res[,1:3],-log(fin_res[,4],base=10),fin_res[,c(5:8,11)],-log(fin_res[,12],base = 10),fin_res[,13:15]))
tot_res[big_to_small_data,-1] <- new_dat
colnames(tot_res) <- c("Big to small ind","Coeff","StdErr","T-value","-log_10(p)","Weight_MMW0.9","WABH_MMW0.9","WABH_MMW0.5","Ten_Rule","CAMT_res","-log_10(p_wabh)","CAMT_pi0","CAMT_k","Weight_MMW0.5")


####################### Creating NIFTI files #####################

Brain_img <- readNIfTI(paste0("Data analysis",delim,"Data",delim,"betsct1_unsmooth.nii"))
tot_res <- data.frame(tot_res)
data_dim <- c(157,189,156)

#### Creating arrays with same dimension as background image:
weig_arr   <- array(0,dim(Brain_img))
camt_pi1   <- array(0,dim(Brain_img))
dis_wabh   <- array(0,dim(Brain_img))
dis_camt   <- array(0,dim(Brain_img))
Tested     <- array(0,dim(Brain_img))
incon_wabh <- array(0,dim(Brain_img))


weig_arr  [1:(157)+12,(1):(189)+14,1:(156)+12]<- array(tot_res$Weight_MMW0.9,data_dim)
camt_pi1  [1:(157)+12,(1):(189)+14,1:(156)+12]<- array(tot_res$CAMT_pi0,data_dim)
dis_wabh  [1:(157)+12,(1):(189)+14,1:(156)+12]<- array(tot_res$WABH_MMW0.9,data_dim)
dis_camt  [1:(157)+12,(1):(189)+14,1:(156)+12]<- array(tot_res$CAMT_res,data_dim)
Tested    [1:(157)+12,(1):(189)+14,1:(156)+12]<- array(1-1*I(tot_res[,4]==0),data_dim)
incon_wabh[1:(157)+12,(1):(189)+14,1:(156)+12]<- array(1*I(tot_res$Weight_MMW0.9<=0.1),data_dim)


weig_arr_nifti  <- oro.nifti::nifti(weig_arr)
camt_pi1_nifti  <- oro.nifti::nifti(camt_pi1)
dis_wabh_nifti  <- oro.nifti::nifti(dis_wabh)
dis_camt_nifti  <- oro.nifti::nifti(dis_camt)
Tested_nifti    <- oro.nifti::nifti(Tested)
incon_wabh_nifti<- oro.nifti::nifti(incon_wabh)

mask = copyNIfTIHeader(img = Brain_img, arr = Tested_nifti)
dis_wabh_nifti = copyNIfTIHeader(img = Brain_img, arr = dis_wabh_nifti)
dis_camt_nifti = copyNIfTIHeader(img = Brain_img, arr = dis_camt_nifti)
incon_wabh_nifti = copyNIfTIHeader(img = Brain_img, arr = incon_wabh_nifti)
weig_arr_nifti = copyNIfTIHeader(img = Brain_img, arr = weig_arr_nifti)
camt_pi1_nifti = copyNIfTIHeader(img = Brain_img, arr = camt_pi1_nifti)

write_nifti(mask, paste0("Data analysis",delim,"Results",delim,"mask.nii"))
write_nifti(dis_wabh_nifti, paste0("Data analysis",delim,"Results",delim,"dis_wabh.nii"))
write_nifti(dis_camt_nifti, paste0("Data analysis",delim,"Results",delim,"dis_camt.nii"))
write_nifti(incon_wabh_nifti, paste0("Data analysis",delim,"Results",delim,"incon_wabh.nii"))
write_nifti(weig_arr_nifti, paste0("Data analysis",delim,"Results",delim,"weig_arr.nii"))
write_nifti(camt_pi1_nifti, paste0("Data analysis",delim,"Results",delim,"camt_pi1.nii"))


