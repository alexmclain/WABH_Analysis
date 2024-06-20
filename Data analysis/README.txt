Title: "False Discovery Rate Control for Lesion-Symptom Mapping with Heterogeneous data via Weighted P-values"

Authors: Siyu Zheng, Alexander C. McLain, Joshua Habiger, Christopher Rorden, and Julius Fridriksson


R version 4.1.1 (2021-08-10)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:

[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252  

[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                         

[5] LC_TIME=English_United States.1252   

attached base packages:

[1] splines   stats     graphics  grDevices utils     datasets  methods   base    

other attached packages:

[1] ggnewscale_0.4.10  neurobase_1.32.3   oro.nifti_0.11.4   adaptMT_1.0.0      data.table_1.14.2

 [6] mgcv_1.8-36        nlme_3.1-152       CAMT_1.1           matrixStats_0.61.0 cowplot_1.1.1    

[11] ggplot2_3.3.5      qvalue_2.26.0      Matrix_1.3-4       dplyr_1.1.1      


loaded via a namespace (and not attached):

[1] Rcpp_1.0.10       pillar_1.9.0      compiler_4.1.1    plyr_1.8.6        R.methodsS3_1.8.1

[6] R.utils_2.11.0    bitops_1.0-7      tools_4.1.1       lifecycle_1.0.3   tibble_3.2.1    

[11] gtable_0.3.0      lattice_0.20-44   pkgconfig_2.0.3   rlang_1.1.0       RNifti_1.6.1    

[16] cli_3.6.1         rstudioapi_0.13   withr_2.4.3       stringr_1.4.0     generics_0.1.3  

[21] vctrs_0.6.1       grid_4.1.1        tidyselect_1.2.0  glue_1.6.2        R6_2.5.1        

[26] fansi_0.5.0       reshape2_1.4.4    magrittr_2.0.1    scales_1.1.1      abind_1.4-5     

[31] colorspace_2.0-2  utf8_1.2.2        stringi_1.7.4     munsell_0.5.0     R.oo_1.24.0 



/folder Data analysis/

01-Logistic-reg.R

An R script which was used to generate the logistic regression model results. The Log_mod_res_logitY_coord_mat.rds 
in the Results subfolder.
Our methods are demonstrated on an aphasia study investigating which regions of the brain are associated with 
the severity of language impairment among stroke survivors. The data for this study are not publically available 
but can be obtained by contacting the study PI's Julius Fridriksson (jfridrik 'at' sc 'dot' edu) or Christopher 
Rorden (RORDEN "at" mailbox 'dot' sc 'dot' edu). For this analysis, we use the results of voxel-by-voxel logistic 
regression analysis of lesion status (yes/no) by Aphasia Quotient (AQ) (logit transformed), and total lesion 
volume (logit transformed).

02a-Adapt.R

An R script which was used to generate part of the results related to the AdaPT estimation methods in Table 1 in the 
Supplemental Material. Adapt_results.rds in the Results subfolder can be used to upload the results of the AdaPT function.

02b-CAMT.R

An R script which was used to generate part of the results related to the CAMT estimation methods in Table 1 in the 
Supplemental Material. CAMT_resuls.rds in the Results subfolder can be used to upload the results of the CAMT function. 
Final_initial_res.rds in the Results subfolder was generated from 02b-CAMT.R Also.

03-Creating_NIfTI_files.R

An R script which use Final_initial_res.rds to generate intermediate results for the data analysis results plotting.

04-Plotting.R

An R script which was used to plot data analysis results on structural brain image.


/subfolder Data/

Y_dat.rds contains all Aphasia Quotient(AQ) data.
betsct1_unsmooth.nii contains all structural brain image data.
big_to_small_data.rds contains indicators which voxels were used in the data analysis.

/subfolder Results/

Log_mod_res_logitY_coord_mat.rds contains all logistic regression model results.
Adapt_results.rds and CAMT_results.rds are the intermediate results for 02a-Adapt.R and 02b-CAMT.R.
Final_initial_res.rds contains all the data analysis results for plotting.
camt_pi1.nii.gz
dis_wabh.nii.gz
incon_wabh.nii.gz
mask.nii.gz
weig_arr.nii.gz are all intermediate results for plotting.
