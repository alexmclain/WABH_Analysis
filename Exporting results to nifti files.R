fin_res <- read.csv("Final_initial_res.csv")
big_to_small_data <- read.csv("Big_to_small_ind.csv")
big_to_small_data <- big_to_small_data[,-1]
sum(1*I(big_to_small_data))

tot_res <- cbind(big_to_small_data,0,0,0,0,0,0,0,0,0,0,0)

tot_res[big_to_small_data,2] <- fin_res[,1]
tot_res[big_to_small_data,3] <- fin_res[,2]
tot_res[big_to_small_data,4] <- fin_res[,3]
tot_res[big_to_small_data,5] <- -log(fin_res[,4],base=10)
tot_res[big_to_small_data,6] <- fin_res[,5]
tot_res[big_to_small_data,7] <- fin_res[,6]
tot_res[big_to_small_data,8] <- fin_res[,7]
tot_res[big_to_small_data,9] <- fin_res[,8]
tot_res[big_to_small_data,10]<- fin_res[,9]
tot_res[big_to_small_data,11]<- -log(fin_res[,11],base = 10)
tot_res[big_to_small_data,12]<- -log(fin_res[,12],base = 10)

colnames(tot_res) <- c("Big to small ind","Coeff","StdErr","T-value","-log_10(p)","Weight","WBonf_res","WABH","Twenty","BH","-log_10(p_wabh)","-log_10(p_10)")
summary(tot_res)


#Creating a 157x189x156 array (correctly).

p_wabh <- array(tot_res[,11],c(157,189,156))
p_20 <- array(tot_res[,12],c(157,189,156))
p_bh <- array(tot_res[,5],c(157,189,156))
weig_arr <- array(tot_res[,6],c(157,189,156))
dis_wabh <- array(tot_res[,8],c(157,189,156))
Tested <- array(1-1*I(tot_res[,4]==0),c(157,189,156))

library("R.matlab")
writeMat("Analysis_output_30Oct2018.mat",p_wabh=p_wabh,p_bh=p_bh,weig_arr=weig_arr,dis_wabh=dis_wabh,Tested=Tested)


#-----------------------------
library(fslr) # need niftiarr
library(oro.nifti)
p_wabh_nifti = oro.nifti::nifti(p_wabh)
p_20_nifti = oro.nifti::nifti(p_20)
p_bh_nifti = oro.nifti::nifti(p_bh)
dis_wabh_nifti = oro.nifti::nifti(dis_wabh)
Incon_dis_wabh_nifti = oro.nifti::nifti(dis_wabh)
weig_arr_nifti = oro.nifti::nifti(weig_arr)
Tested_nifti = oro.nifti::nifti(Tested)


mask = niftiarr(Tested_nifti,Tested_nifti==1)
p_wabh_nifti[mask == 0] = 0
p_bh_nifti[mask == 0] = 0
p_20_nifti[mask == 0] = 0
dis_wabh_nifti[mask == 0] = 0
Incon_dis_wabh_nifti[-weig_arr_nifti>=(-0.1)] = 2
Incon_dis_wabh_nifti[mask == 0] = 0
weig_arr_nifti[mask == 0] = 0
mask[mask == 0] = 0





####################### Exporting NIFTI files #####################


p_wabh <- array(tot_res[,11],c(157,189,156))
p_20 <- array(tot_res[,12],c(157,189,156))
p_bh <- array(tot_res[,5],c(157,189,156))
weig_arr <- array(tot_res[,6],c(157,189,156))
dis_wabh <- array(tot_res[,8],c(157,189,156))
#Tested <- array(1-1*I(tot_res[,4]==0),c(157,189,156))

rows_out <- 25

p_wabh[,,1:(dim(p_wabh)[3]-rows_out)] <- p_wabh[,,(rows_out+1):(dim(p_wabh)[3])]
p_wabh[,,(dim(p_wabh)[3]-rows_out+1):(dim(p_wabh)[3])] <- 0
p_20[,,1:(dim(p_wabh)[3]-rows_out)] <- p_20[,,(rows_out+1):(dim(p_wabh)[3])]
p_20[,,(dim(p_wabh)[3]-rows_out+1):(dim(p_wabh)[3])] <- 0
p_bh[,,1:(dim(p_wabh)[3]-rows_out)] <- p_bh[,,(rows_out+1):(dim(p_wabh)[3])]
p_bh[,,(dim(p_wabh)[3]-rows_out+1):(dim(p_wabh)[3])] <- 0
dis_wabh[,,1:(dim(p_wabh)[3]-rows_out)] <- dis_wabh[,,(rows_out+1):(dim(p_wabh)[3])]
dis_wabh[,,(dim(p_wabh)[3]-rows_out+1):(dim(p_wabh)[3])] <- 0
weig_arr[,,1:(dim(p_wabh)[3]-rows_out)] <- weig_arr[,,(rows_out+1):(dim(p_wabh)[3])]
weig_arr[,,(dim(p_wabh)[3]-rows_out+1):(dim(p_wabh)[3])] <- 0
#Tested[,,1:(dim(p_wabh)[3]-rows_out)] <- Tested[,,(rows_out+1):(dim(p_wabh)[3])]
#Tested[,,(dim(p_wabh)[3]-rows_out+1):(dim(p_wabh)[3])] <- 0



p_wabh_nifti = oro.nifti::nifti(p_wabh)
p_20_nifti = oro.nifti::nifti(p_20)
p_bh_nifti = oro.nifti::nifti(p_bh)
weig_arr_nifti = oro.nifti::nifti(weig_arr)
dis_wabh_nifti = oro.nifti::nifti(dis_wabh)
#Tested_nifti = oro.nifti::nifti(Tested)
Incon_dis_wabh_nifti = oro.nifti::nifti(dis_wabh*0)




#mask = niftiarr(Tested_nifti,Tested_nifti==1)
p_wabh_nifti[mask == 0] = -1
p_bh_nifti[mask == 0] = -1
p_20_nifti[mask == 0] = -1
dis_wabh_nifti[mask == 0] = 0
Incon_dis_wabh_nifti[-weig_arr_nifti>=(-0.1)] = 2
Incon_dis_wabh_nifti[mask == 0] = 0
weig_arr_nifti[mask == 0] = -1
#mask[mask == 0] = 0

p_wabh_nifti = p_wabh_nifti+1
p_bh_nifti = p_bh_nifti+1
weig_arr_nifti= weig_arr_nifti +1

sWeights <- readNIfTI("Final_results/alex/sWeights.nii")
#mask = copyNIfTIHeader(img = sWeights, arr = mask)
p_wabh_nifti = copyNIfTIHeader(img = sWeights, arr = p_wabh_nifti)
p_bh_nifti = copyNIfTIHeader(img = sWeights, arr = p_bh_nifti)
dis_wabh_nifti = copyNIfTIHeader(img = sWeights, arr = dis_wabh_nifti)
Incon_dis_wabh_nifti = copyNIfTIHeader(img = sWeights, arr = Incon_dis_wabh_nifti)
weig_arr_nifti = copyNIfTIHeader(img = sWeights, arr = weig_arr_nifti)

show(weig_arr_nifti)

#writeNIfTI2(mask,"Final_results/Mask.nii.gz")
writeNIfTI2(p_wabh_nifti,"Final_results/P_values_weight.nii.gz")
writeNIfTI2(p_bh_nifti,"Final_results/P_values_std.nii.gz")
writeNIfTI2(dis_wabh_nifti,"Final_results/Discoveries.nii.gz")
writeNIfTI2(Incon_dis_wabh_nifti,"Final_results/DiscIncon.nii.gz")
writeNIfTI2(weig_arr_nifti,"Final_results/Weights.nii.gz")


# The dimension of each figure should be c(157,189,156), which were the dimensions of the original data.  Here is a discription of file is each:
#    - Mask: this is a 1 if the voxel produced a p-value, NA otherwise. A voxel produced a p-value if there was at least 1 lesioned patient and the logistic algorithm converged.
#    - P_values_std: -log_10(P_values).  No adjustement made.
#    - P_values_weight: -log_10(Q_values), where the weighting procedure was applied. 
#    - Discoveries: 1 if the voxel was found to be significant (using the weighted BH procedure), 0 otherwise.
#    - Weights: the p-values weights

par(mfrow=c(5,5),mai=c(0,0,0,0))

for(k in 83:106){
d <- dim(mask)
mat <- expand.grid(1:c(d[1]),1:c(d[2]))
plot(mat[,1],mat[,2],col=1,pch=15,axes=FALSE,ann=FALSE,xlim = c(1,d[1]),ylim = c(1,d[2]))

vals <- round(array(mask[,,k]),3)
col_v <- rep(1,length(vals))

t1 <- vals
col_v[vals==1] <- 3

points(mat[,1],mat[,2],col=col_v,pch=15)
}


rect <- locator()

for(j in 1:length(rect$x)){mask[floor(rect$x[j]):d[1],floor(rect$y[j]):d[2],k] <- 0}


par(mfrow=c(1,1),mai=c(0,0,0,0))

for(k in 83:106){
  d <- dim(p_bh_nifti)
  mat <- expand.grid(1:c(d[1]),1:c(d[2]))
  plot(mat[,1],mat[,2],col=1,pch=15,axes=FALSE,ann=FALSE,xlim = c(1,d[1]),ylim = c(1,d[2]))
  
  vals <- array(Incon_dis_wabh_nifti[,,k])
  col_v <- rep(1,length(vals))
  
  t1 <- vals
  col_v <- gray(vals/max(vals))
  
  points(mat[,1],mat[,2],col=col_v,pch=15)
}


rect <- locator()

for(j in 1:length(rect$x)){mask[floor(rect$x[j]):d[1],floor(rect$y[j]):d[2],k] <- 0}

