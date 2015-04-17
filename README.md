---
title: "F2_EW_MainAnalysis_README"
author: "Quentin Schorpp"
date: "2015-04-17"
output: md_document
---

*The README.md file gives an overview of all the files in the project. It should briefly describe the project including things like its title, author(s), topic, any copyright information, and so on. It should also indicate how the folders in the project are organized and give instruction how to reproduce the project. The README file should be in the main project folder. It is good practice to dynamically include the system information for the R session that was used. To do this you can write your README file with R Markdown. Simply include the sessionInfo() command in a knitr code chunk in the R Markdown document. If you knit this file immediately after knitting your presentation document it will record the information for that session.*

## Topic: 
Earthworm communities from fields of the perennial bioenergy crop *S. perfoliatum*.

## Project Directories:
**Data:** Main Data Source (F2_EW_Data.xlsx; F2_EW_Total.txt)
      Subsets of Data from F2_EW_Total.txt (other .txt files)
      Files: 
      
      * F2_EW_Total.txt - Dataset with all variables and all parameters
      * F2_EW_Data.xlsx - Excel spreadsheet for porucing .txt files and subsets
      
**Data/GatherSource:** Data Processing 
                    Files:
                    
                    * Gather1.R - Summing up for functional groups; 
                                  Calculation of Biodiversity Indices;
                                  1rst Order Data [180 rows]
                                  
                    
**Analysis:** Main Analysis
           Files:
           
           *F2_EW_MainAnalysis.Rmd - Analysis of earthworm abundances and biomass using glmm




```
## R version 3.1.3 (2015-03-09)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 7 x64 (build 7601) Service Pack 1
## 
## locale:
## [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
## [3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
## [5] LC_TIME=German_Germany.1252    
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] extrafont_0.17    scales_0.2.4      plyr_1.8.1       
##  [4] reshape2_1.4.1    coefplot2_0.1.3.2 coda_0.17-1      
##  [7] lme4_1.1-7        Rcpp_0.11.5       Matrix_1.2-0     
## [10] bbmle_1.0.17      glmmADMB_0.7.7    R2admb_0.7.13    
## [13] MASS_7.3-40       ggplot2_1.0.1    
## 
## loaded via a namespace (and not attached):
##  [1] cluster_2.0.1     colorspace_1.2-6  digest_0.6.8     
##  [4] evaluate_0.6      extrafontdb_1.0   formatR_1.1      
##  [7] gtable_0.1.2      htmltools_0.2.6   knitr_1.9        
## [10] labeling_0.3      lattice_0.20-31   mgcv_1.8-6       
## [13] minqa_1.2.4       munsell_0.4.2     nlme_3.1-120     
## [16] nloptr_1.0.4      numDeriv_2012.9-1 parallel_3.1.3   
## [19] permute_0.8-3     proto_0.3-10      reshape_0.8.5    
## [22] rmarkdown_0.5.1   Rttf2pt1_1.3.3    splines_3.1.3    
## [25] stringr_0.6.2     tools_3.1.3       vegan_2.2-1      
## [28] yaml_2.1.13
```
