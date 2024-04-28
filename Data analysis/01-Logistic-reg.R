### Our methods are demonstrated on an aphasia study investigating which regions of 
### the brain are associated with the severity of language impairment among stroke 
### survivors. The data for this study are not publically available but can be obtained 
### by contacting the study PI's Julius Fridriksson (jfridrik 'at' sc 'dot' edu) or 
### Christopher Rorden (RORDEN "at" mailbox 'dot' sc 'dot' edu). For this analysis, 
### we use the results of voxel-by-voxel logistic regression analysis of lesion status 
### (yes/no) by Aphasia Quotient (AQ) (logit transformed), and total lesion volume 
### (logit transformed). Please see the full manuscript for more details.

## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

library(dplyr)
library(Matrix)
library(qvalue)


CD20    <- readRDS("Data/Individual_data.rds") # Individual level data
Big_img <- readRDS("Data/Big_img.rds") # Full image data
coord_mat <- readRDS("Data/coord_mat.rds") # 3d Coordinates of image data
zero_vec <- readRDS("Data/zero_vec.rds") # Locations with zero damage

## Getting valid rows
valid_rows <- c(1:length(zero_vec))[as.vector(zero_vec)==0]
M <- length(valid_rows)
coord_mat <- coord_mat[valid_rows,]

remove(zero_vec)

Y_dat <- CD20$Y_dat #Outcome WAB Score
X_pl <- CD20$Xp #Total lesion size
X_logit <- log(X_pl/M/(1-X_pl/M)) #Logis transform

pdf("Figures/Xim_data_analysis.pdf", width = 5, height = 5)
hist(X_logit, breaks = 40, 
     main = expression(Data~Analysis~X[im]^'+'),
     xlab = expression(X[im]^'+'))
dev.off()

# Reverse logit transform of Y 
Y_dat <- (100-Y_dat)/100
Y_dat <- log(Y_dat)
N <- length(Y_dat)

# Model fit
plot(Y_dat,X_logit)


#Run the analysis and save the results (might take a while).
X_mean <- rep(0,M)
R_sq <- rep(0,M)
lm_info <- matrix(0,M,4)
for(i in 1:M){
  BR_dat <- c(as.vector(Big_img[valid_rows[i],]))
  BR_dat <- BR_dat[complete.cases(Cov_dat)]
  X_mean[i] <- mean(BR_dat)
  if(X_mean[i]==0){
    t_lm_info <- c(0,1,0,1)
    cat(i,"zero mean error \n")}
  if(X_mean[i]>0){
    t_X_pl <- (X_pl-BR_dat)/(M-1)
    t_X_pl <- log(t_X_pl/(1-t_X_pl))
    t_X_pl <- (t_X_pl - mean(t_X_pl))/sd(t_X_pl)
    t1 <- try(t_lm <- lm(Y_dat~t_X_pl))
    if(attr(t1,"class")[1] == "try-error"){
      t_lm_info <- c(0,1,0,1)
      cat(i,"LM error \n")}
    if(attr(t1,"class")[1] == "lm"){
      R_sq[i] <- summary(t_lm)$r.squared
      t1 <- try(t_glm <- glm(BR_dat ~Y_dat+t_X_pl,family=binomial(link='logit')),silent = TRUE)
      if(attr(t1,"class")[1] == "try-error"){
        t_lm_info <- c(0,1,0,1)
        cat(i,"GLM error \n")}
      if(attr(t1,"class")[1] == "glm"){
        t_lm_info <- c(summary(t_glm)$coefficients[2,1:4])}
    }
  }
  lm_info[i,] <- t_lm_info
}

# Transform to one-tailed p-values
lm_info[,4] <- pnorm(lm_info[,3],lower.tail = FALSE)

# Save results
res <- data.frame(X_mean,R_sq,lm_info,coord_mat)
saveRDS(res,paste0("Data analysis",delim,"Results",delim,"Log_mod_res_logitY_coord_mat.rds"))

