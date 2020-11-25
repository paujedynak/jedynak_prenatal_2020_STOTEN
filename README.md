# Prenatal exposure to a wide range of environmental chemicals and child behaviour between 3 and 7 years of age - An exposome-based approach in 5 European cohorts

Repository contains reproducible analyses for the paper: Jedynak et al., 2020 (**DOI: https://doi.org/10.1016/j.scitotenv.2020.144115**). Depends on R (>= 3.5.0).


## Analysis overview

This analysis was performed under Windows 10 x64 (build 19041) using:    
* [R 4.0.2](https://cran.r-project.org/bin/windows/base) (2020-06-22)   
* [RStudio 1.3.1056](https://rstudio.com)    
* [renv](https://cran.r-project.org/web/packages/renv/index.html) dependency management package    
* [drake](https://github.com/ropensci/drake) workflow management package    


### To run

Re-running the analysis requires executing of `renv::hydrate()` or `renv::restore()` which will upgrade/ downgrade user's packages to the versions used in the present study. This operation will modify the packages only locally (for the project), so it will not affect user's package library.

Re-running the analysis requires an additional `data/raw_data` folder not shared here, containing the different data sets presented hereafter. These data can only be provided upon request and after approval by the Helix consortium.

```bash
Rscript makefile.R
```

Only the `makefile.R` script needs to be executed. This makefile calls the drake analysis plan (`plan.R`) which outlines the successive analysis steps and calls the different functions used to execute different analysis steps, to be found in the `R/` folder. 


## Repo organization

Analysis input data-files are not made available as they contain sensitive information. The analysis input data files are:

* `SDQ_biomarkers_july_2018.csv` = dataset containing exposures and SDQ scores data
* `20180212_v2_3_ImputedDataset_20x10_AllDataFrame.RDS` = dataset containing covariates data
* `INMA_EDEN_KANK_RHEA_BIB_SDQ_2018_07.csv` = dataset containing sex and age of the children


### data/ folder

Analysis input data-files are not made available as they contain sensitive information. Files in the `results/variable_lists` folder contain manually created lists of variables.


### results/ folder

Results files in the `results/` folder:

* `results/main_analyses` folder contains results that were used in the paper in the manuscript.
* `results/main_analyses/tables` folder contains tables from the paper in .csv format.
* `results/main_analyses/figures` folder contains figures from the paper in .jpg format.    

* `results/supplementary_analyses` folder contains results that were supplementary or mentioned in the paper as "not shown" - these are sensitivity analyses for additional adjustment of the ExWAS for breastfeeding and for fish and seafood consumption.    

* `results/revision` folder contains additional analyses that were included in the answer to the Reviewer.

Code used to produce the files is found in the corresponding .Rmd files (jedynak_prenatal_2020.Rmd for all analyses for the paper and jedynak_prenatal_2020_answer_to_reviewer.Rmd for the answer to the Reviewer).


### R/ folder

This folder contains all the packages and functions used for the analyses, called in the `makefile.R` script.

* `packages.R` = analysis packages.
* `variables_lists.R` = miscellaneous helper lists of compounds, variables, etc.
* `plan.R` = drake plan outlying the different steps of the analysis and setting the dependencies.
* `tables.R` = code for all paper tables
* `figures.R` = code for all paper figures.


### Packages

All packages used in the analyses are saved in the `renv/library/R-4.0/x86_64-w64-mingw32` folder at the version used to produce these results, under control of the `renv` package manager.


## Session info

```
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] PrepareData_0.1.0        MultipleRegression_0.1.0 MultipleImputation_0.1.0
 [4] LassoENET_0.1.0          Heterogeneity_0.1.0      Helpers_0.1.0           
 [7] GAMrcs_0.1.0             ExWAS_0.1.0              DescriptiveStats_0.1.0  
[10] forcats_0.5.0            stringr_1.4.0            purrr_0.3.4             
[13] tidyr_1.1.2              tibble_3.0.4             tidyverse_1.3.0         
[16] readr_1.4.0              mice_3.11.0              magrittr_1.5            
[19] MASS_7.3-53              lubridate_1.7.9          lme4_1.1-23             
[22] Matrix_1.2-18            Hmisc_4.4-1              ggplot2_3.3.2           
[25] Formula_1.2-3            survival_3.2-7           lattice_0.20-41         
[28] here_0.1                 haven_2.3.1              drake_7.12.6            
[31] dplyr_1.0.2              broom_0.7.1             

loaded via a namespace (and not attached):
  [1] rms_6.0-1           tidyselect_1.1.0    htmlwidgets_1.5.2   grid_4.0.3         
  [5] munsell_0.5.0       base64url_1.4       codetools_0.2-16    chron_2.3-56       
  [9] statmod_1.4.34      withr_2.3.0         colorspace_1.4-1    filelock_1.0.2     
 [13] knitr_1.30          uuid_0.1-4          rstudioapi_0.11     stats4_4.0.3       
 [17] pscl_1.5.5          officer_0.3.14      bbmle_1.0.23.1      txtq_0.2.3         
 [21] rprojroot_1.3-2     vctrs_0.3.4         generics_0.0.2      TH.data_1.0-10     
 [25] metafor_2.4-0       xfun_0.18           R6_2.4.1            doParallel_1.0.15  
 [29] assertthat_0.2.1    scales_1.1.1        multcomp_1.4-14     nnet_7.3-14        
 [33] gtable_0.3.0        conquer_1.0.2       sandwich_3.0-0      rlang_0.4.8        
 [37] MatrixModels_0.4-1  systemfonts_0.3.2   Rmisc_1.5           splines_4.0.3      
 [41] checkmate_2.0.0     yaml_2.2.1          abind_1.4-5         modelr_0.1.8       
 [45] backports_1.1.10    HardyWeinberg_1.6.8 tools_4.0.3         ellipsis_0.3.1     
 [49] kableExtra_1.2.1    RColorBrewer_1.1-2  Rsolnp_1.16         Rcpp_1.0.5         
 [53] plyr_1.8.6          visNetwork_2.0.9    base64enc_0.1-3     progress_1.2.2     
 [57] prettyunits_1.1.1   rpart_4.1-15        zoo_1.8-8           cluster_2.1.0      
 [61] fs_1.5.0            data.table_1.13.0   forestplot_1.10     openxlsx_4.2.2     
 [65] flextable_0.5.11    SparseM_1.78        reprex_0.3.0        truncnorm_1.0-8    
 [69] mvtnorm_1.1-1       storr_1.2.4         matrixStats_0.57.0  hms_0.5.3          
 [73] evaluate_0.14       XML_3.99-0.5        rio_0.5.16          jpeg_0.1-8.1       
 [77] readxl_1.3.1        gridExtra_2.3       Gmisc_1.11.0        MAc_1.1.1          
 [81] compiler_4.0.3      bdsmatrix_1.3-4     writexl_1.3.1       crayon_1.3.4       
 [85] minqa_1.2.4         htmltools_0.5.0     mpath_0.3-26        DBI_1.1.0          
 [89] dbplyr_1.4.4        boot_1.3-25         car_3.0-10          cli_2.1.0          
 [93] parallel_4.0.3      igraph_1.2.6        pkgconfig_2.0.3     metaplus_0.7-11    
 [97] numDeriv_2016.8-1.1 foreign_0.8-80      xml2_1.3.2          foreach_1.5.0      
[101] webshot_0.5.2       rvest_0.3.6         digest_0.6.25       fastGHQuad_1.0     
[105] rmarkdown_2.4       cellranger_1.1.0    gam_1.20            htmlTable_2.1.0    
[109] gdtools_0.2.2       curl_4.3            quantreg_5.73       compareGroups_4.4.5
[113] nloptr_1.2.2.2      lifecycle_0.2.0     nlme_3.1-149        jsonlite_1.7.1     
[117] carData_3.0-4       viridisLite_0.3.0   fansi_0.4.1         pillar_1.4.6       
[121] bst_0.3-21          httr_1.4.2          glue_1.4.2          zip_2.1.1          
[125] gbm_2.1.8           png_0.1-7           iterators_1.0.12    stringi_1.5.3      
[129] blob_1.2.1          polspline_1.1.19    latticeExtra_0.6-29 renv_0.12.0    

```

