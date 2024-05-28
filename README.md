# WABH_Analysis
This project considers the use of p-value weighting to voxel-based lesion behavior mapping (VLSM) studies. The programs to estimate the MMW weights proposed in the paper are available in the **Functions** folder.  The estimation algorithm requires an estimate of the non-null probability. This can be done by test or overall. See the analysis scripts below for a demonstration of the methods.

Our methods are demonstrated on an aphasia study investigating which regions of the brain are associated with the severity of language impairment among stroke survivors. The data for this study are not publically available but can be obtained by contacting the study PI's Julius Fridriksson (jfridrik 'at' sc 'dot' edu) or Christopher Rorden (RORDEN "at" mailbox 'dot' sc 'dot' edu). 

See the **Data Analysis** folder (and the corresponding **Read me**), which contains all the `R` scripts to reproduce the data example in the main paper. For this analysis, we use the results of voxel-by-voxel logistic regression analysis of lesion status (yes/no) by Aphasia Quotient (AQ) (logit transformed), and total lesion volume (logit transformed). While the full data are not publically available, we do provide the voxel-by-voxel results of the logistic regression analyses. From these results, all multiple testing methods and results can be applied.

The **Simulation Studies** folder has all the functions to reproduce the iteration-by-iteration results, along with the iteration-by-iteration and summarized results. Code to reproduce all figures and tables for the simulations are included.






