# WABH_Analysis
This project considers the use of p-value weighting to voxel-based lesion behavior mapping (VLSM) studies. The programs to estimate the MMW weights proposed in the paper are available in the **Programs** folder.  The estimation algorithm requires an estimate of the non-null probability. This can be done by test or overall. See the analysis scripts below for a demonstration of the methods.

Our methods are demonstrated on an aphasia study investigating which regions of the brain are associated with the severity of language impairment among stroke survivors. The data for this study are not publically available but can be obtained by contacting the study PI's Julius Fridriksson (jfridrik 'at' sc 'dot' edu) or Christopher Rorden (RORDEN "at" mailbox 'dot' sc 'dot' edu). For this analysis, we use the results of voxel-by-voxel logistic regression analysis of lesion status (yes/no) by Aphasia Quotient (AQ) (logit transformed), and total lesion volume (logit transformed). Please see the full manuscript for more details. 

The data are available in the **Data** folder, the scripts to reproduce the results in the paper are available in the **Analysis scripts** folder.


