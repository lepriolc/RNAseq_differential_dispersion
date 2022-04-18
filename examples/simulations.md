# Examples of parameters for the analysis of simulated RNA-seq datasets

This file contains examples of parameters and R session informations for the generation and the analysis of simulated RNA-seq datasets.

## Generate datasets

This section contains an example of parameter section to use in the analysis script `simulations-generateDatasets.R` and the R session information to generate a simulated RNA-seq dataset.

### Parameters

```
##############
# Parameters #
##############
# number of samples per condition
samples_per_condition <- 50
# number of genes
genes <- 10000
# sequencing depth
depth <- 5000000
min_depth_param <- max_depth_param <- 1
# DE genes
DE_fraction <- 0
# DE_fraction <- 0
DE_down_fraction <-  0.5 # fraction of DE with a decrease of mean in the second condition
DE_base_effect <- 1.5
non_DE_runif_FC <- TRUE # generate FC for non DE genes using runif() function
# DD genes
DD_fraction <- 0.5
DD_down_fraction <- 0.5 # fraction of DD with a decrease of dispersion in the second condition
DD_base_effect <- 1.5
DD_effect_param <- 1
non_DD_runif_FC <- FALSE # generate FC for non DD genes using runif() function
# outliers
single_outlier_high_fraction <- 0.1
single_outlier_low_fraction <- 0
random_outlier_high_fraction <- 0
random_outlier_low_fraction <- 0
# number of replicates
replicates <- 10
# path to output directory
output_dir <- "./output/simulations/00-Data"
```

This setting generates 10 replicates of a dataset composed of 2 sets of 50 samples. It contains only lowly differentially expressed genes with a miximum of 1.5 fold-change between the sets of samples.

The datasets analyzed in the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression' were generated with the following parameter values:

- samples_per_condition = {20; 30; 40; 50; 100}
- DE_fraction = {0; 0.5}
- DE_base_effect = {1.1; 1.2; 1.3; 1.4; 1.5}

### R session information

```
R version 3.6.2 (2019-12-12)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/openblas-base/libblas.so.3
LAPACK: /usr/lib/libopenblasp-r0.2.18.so

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.3.4    compcodeR_1.20.1 sm_2.2-5.6      

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3          compiler_3.6.2      pillar_1.4.2        plyr_1.8.4          timeSeries_3062.100
 [6] bitops_1.0-6        tools_3.6.2         vioplot_0.3.4       rpart_4.1-15        clue_0.3-57        
[11] stable_1.1.4        tibble_2.1.3        gtable_0.2.0        lattice_0.20-38     pkgconfig_2.0.1    
[16] rlang_0.4.11        rstudioapi_0.7      rmutil_1.1.3        xfun_0.29           statip_0.2.3       
[21] withr_2.1.2         stringr_1.3.1       cluster_2.1.0       dplyr_0.8.3         knitr_1.37         
[26] gtools_3.8.1        caTools_1.17.1      locfit_1.5-9.1      grid_3.6.2          tidyselect_0.2.5   
[31] modeest_2.4.0       glue_1.3.1          R6_2.2.2            tcltk_3.6.2         spatial_7.3-11     
[36] gdata_2.18.0        limma_3.36.2        purrr_0.2.5         ROCR_1.0-7          edgeR_3.22.3       
[41] magrittr_1.5        scales_0.5.0        gplots_3.0.1        fBasics_3042.89     MASS_7.3-50        
[46] stabledist_0.7-1    assertthat_0.2.0    timeDate_3043.102   colorspace_1.3-2    KernSmooth_2.23-15 
[51] stringi_1.2.3       munsell_0.5.0       markdown_0.8        crayon_1.3.4        zoo_1.8-6
```


## Differential dispersion analysis

This section contains an example of parameter section to use in the analysis script `simulations-DD_analysis.R` and the R session information to identify differentially dispersed genes in a simulated RNA-seq dataset. The parameter setting enables the anlalysis of a replicate generated with the parameter setting of the previous section.

### Parameters

```
##############
# Parameters #
##############
# simulated dataset
## rds file
rds_file <- "./output/simulations/00-Data/10-lowly_DE/base_effect_1_5/50/repl_1/10000_50_50_5000000_DE_0_1_5_0_runif_DD_0_5_1_5_1_1_SO_0_1_0_RO_0_0_repl_1.rds"
## true dispersion annotation column basename
dispersion_colname <- "truedispersions"
## true mean annotation column basename
mean_colname <- "truemeans"
# normalization method
normalization_method <- "TMM"
# minimum expression threshold in CPM to filter out lowly expressed genes
filter_threshold <- 1
# MDSeq parameters
## perform MDSeq outlier removal function
outlier_removal <- TRUE
## fold-change threshold to identify DD genes (not in log2 scale)
DD_FC_threshold <- 1
## fold-change threshold to identify highly genes (not in log2 scale)
DE_FC_threshold <- 1.5
## p-value threshold
pval_threshold <- 0.05
## number of cores to use for computations
cores <- 4
## minimum sample size
min_sample_size <- 1
# path to output directory
# output_dir <- "./output"
output_dir <- "./output/simulations"
```

The path to each replicate of simulated datasets must be provided separately using the 'rds_file' variable.

The datasets analyzed in the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression' were generated with the following parameter values:

- DE_FC_threshold = {1.1; 1.2; 1.3; 1.4; 1.5}

### R session information

```
R version 3.6.2 (2019-12-12)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/openblas-base/libblas.so.3
LAPACK: /usr/lib/libopenblasp-r0.2.18.so

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ROCR_1.0-7     gplots_3.0.1   MDSeq_1.0.5    DiPhiSeq_0.2.0 edgeR_3.22.3   limma_3.36.2  

loaded via a namespace (and not attached):
 [1] compiler_3.6.2     tools_3.6.2        Rcpp_1.0.3         KernSmooth_2.23-15 gdata_2.18.0      
 [6] grid_3.6.2         locfit_1.5-9.1     caTools_1.17.1     bitops_1.0-6       gtools_3.8.1      
[11] lattice_0.20-38
```

