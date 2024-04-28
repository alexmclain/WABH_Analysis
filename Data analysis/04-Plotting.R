## Begin by setting the working drive to location to where the
## folders "Data Analysis", "Functions", and "Simulation Study" are.

delim <- "/"

library(neurobase)
library(oro.nifti)
library(ggplot2)
library(ggnewscale)

####################### Read in NIFTI files #####################

Brain_img <- readNIfTI(paste0("Data analysis",delim,"Data",delim,"betsct1_unsmooth.nii"))
mask <- readNIfTI(paste0("Data analysis",delim,"Results",delim,"mask.nii"))


# Figure of discoveries for WABH method
dis_wabh_nifti <- readNIfTI(paste0("Data analysis",delim,"Results",delim,"dis_wabh.nii"))

z_seq = floor(seq(58,114,length.out = 5))
# Extract matching axial slices
slice_num = z_seq[1]  
anatomical_slice = Brain_img[,,slice_num]
functional_slice = dis_wabh_nifti[,,slice_num]
# Convert matrix slices to data frames for ggplot
anat_df <- expand.grid(x = 1:nrow(anatomical_slice), y = 1:ncol(anatomical_slice), slice = slice_num)
anat_df$intensity <- as.vector(anatomical_slice)
anat_df$signal <- as.vector(functional_slice)

# loop over remaining slices
for(slice_num in z_seq[-1]){
  anatomical_slice = Brain_img[,,slice_num]
  functional_slice = dis_wabh_nifti[,,slice_num]

  t_anat_df <- anat_df  
  t_anat_df$slice <- slice_num
  t_anat_df$intensity <- as.vector(anatomical_slice)
  t_anat_df$signal <- as.vector(functional_slice)
  
  anat_df <- rbind(anat_df, 
                   t_anat_df)
} 

anat_df$x <- 182-anat_df$x
anat_df <- anat_df[anat_df$x<91,]

# Create the base anatomical plot
p <- ggplot() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient(low = "black", high = "white") + 
  new_scale_fill() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = signal)) +
  scale_fill_gradient(low = "transparent", high = "red") +
  coord_fixed() +
  theme_minimal() +
  facet_grid(~slice) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        title = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),  
        panel.background = element_blank(),
        legend.position = "none", 
        panel.spacing = unit(0, "lines"))
p
ggsave(paste0("Data analysis",delim,"Figures",delim,"dis_wabh.pdf"),plot = p, width = 6, height = 3)




# Figure of inconclusive for WABH method
incon_wabh_nifti <- readNIfTI(paste0("Data analysis",delim,"Results",delim,"incon_wabh.nii"))
# voxels not tested are not plotted as inconclusive.
incon_wabh_nifti[mask==0] <- 0

z_seq = floor(seq(58,114,length.out = 5))
# Extract matching axial slices
slice_num = z_seq[1]  
anatomical_slice = Brain_img[,,slice_num]
functional_slice = incon_wabh_nifti[,,slice_num]
# Convert matrix slices to data frames for ggplot
anat_df <- expand.grid(x = 1:nrow(anatomical_slice), y = 1:ncol(anatomical_slice), slice = slice_num)
anat_df$intensity <- as.vector(anatomical_slice)
anat_df$signal <- as.vector(functional_slice)

# loop over remaining slices
for(slice_num in z_seq[-1]){
  anatomical_slice = Brain_img[,,slice_num]
  functional_slice = incon_wabh_nifti[,,slice_num]
  
  t_anat_df <- anat_df  
  t_anat_df$slice <- slice_num
  t_anat_df$intensity <- as.vector(anatomical_slice)
  t_anat_df$signal <- as.vector(functional_slice)
  
  anat_df <- rbind(anat_df, 
                   t_anat_df)
} 

anat_df$x <- 182-anat_df$x
anat_df <- anat_df[anat_df$x<91,]

# Create the base anatomical plot
p <- ggplot() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient(low = "black", high = "white") + 
  new_scale_fill() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = signal)) +
  scale_fill_gradient(low = "transparent", high = "blue") +
  coord_fixed() +
  theme_minimal() +
  facet_grid(~slice) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        title = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),  
        panel.background = element_blank(),
        legend.position = "none", 
        panel.spacing = unit(0, "lines"))
p
ggsave(paste0("Data analysis",delim,"Figures",delim,"incon_wabh.pdf"),plot = p, width = 6, height = 3)





# Figure of weights for WABH method
weig_arr_nifti <- readNIfTI(paste0("Data analysis",delim,"Results",delim,"weig_arr.nii"))
weig_arr_nifti[mask==0] <- NA

z_seq = floor(seq(58,114,length.out = 5))
# Extract matching axial slices
slice_num = z_seq[1]  
anatomical_slice = Brain_img[,,slice_num]
functional_slice = weig_arr_nifti[,,slice_num]
# Convert matrix slices to data frames for ggplot
anat_df <- expand.grid(x = 1:nrow(anatomical_slice), y = 1:ncol(anatomical_slice), slice = slice_num)
anat_df$intensity <- as.vector(anatomical_slice)
anat_df$signal <- as.vector(functional_slice)

# loop over remaining slices
for(slice_num in z_seq[-1]){
  anatomical_slice = Brain_img[,,slice_num]
  functional_slice = weig_arr_nifti[,,slice_num]
  
  t_anat_df <- anat_df  
  t_anat_df$slice <- slice_num
  t_anat_df$intensity <- as.vector(anatomical_slice)
  t_anat_df$signal <- as.vector(functional_slice)
  
  anat_df <- rbind(anat_df, 
                   t_anat_df)
} 

anat_df$x <- 182-anat_df$x
anat_df <- anat_df[anat_df$x<91,]


ggplot() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = signal)) +
  scale_fill_gradient(na.value = "transparent", low = "purple", high = "yellow") +
  coord_fixed() +
  theme_minimal() +
  facet_grid(~slice) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        title = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),  
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"))


# Create the base anatomical plot
p <- ggplot() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") + 
  new_scale_fill() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = signal)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  coord_fixed() +
  theme_minimal() +
  facet_grid(~slice) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        title = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),  
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"))
p
ggsave(paste0("Data analysis",delim,"Figures",delim,"weight_wabh.pdf"),plot = p, width = 6, height = 2.7)





# Figure of pi_hat for the CAMT method
camt_pi1_nifti <- readNIfTI(paste0("Data analysis",delim,"Results",delim,"camt_pi1.nii"))
camt_pi1_nifti[mask==0] <- NA

z_seq = floor(seq(58,114,length.out = 5))
# Extract matching axial slices
slice_num = z_seq[1]  
anatomical_slice = Brain_img[,,slice_num]
functional_slice = camt_pi1_nifti[,,slice_num]
# Convert matrix slices to data frames for ggplot
anat_df <- expand.grid(x = 1:nrow(anatomical_slice), y = 1:ncol(anatomical_slice), slice = slice_num)
anat_df$intensity <- as.vector(anatomical_slice)
anat_df$signal <- as.vector(functional_slice)

# loop over remaining slices
for(slice_num in z_seq[-1]){
  anatomical_slice = Brain_img[,,slice_num]
  functional_slice = camt_pi1_nifti[,,slice_num]
  
  t_anat_df <- anat_df  
  t_anat_df$slice <- slice_num
  t_anat_df$intensity <- as.vector(anatomical_slice)
  t_anat_df$signal <- as.vector(functional_slice)
  
  anat_df <- rbind(anat_df, 
                   t_anat_df)
} 

anat_df$x <- 182-anat_df$x
anat_df <- anat_df[anat_df$x<91,]

# Create the base anatomical plot
p <- ggplot() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") + 
  new_scale_fill() +
  geom_tile(data = anat_df, aes(x = x, y = y, fill = signal)) +
  scale_fill_distiller(palette = "YlGn", na.value = "transparent", direction = -1) +
  coord_fixed() +
  theme_minimal() +
  facet_grid(~slice) +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        title = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),  
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines")
        )
p
ggsave(paste0("Data analysis",delim,"Figures",delim,"pihat_camt.pdf"),plot = p, width = 6, height = 2.7)

