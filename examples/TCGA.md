# Examples of parameters for the analysis of TCGA RNA-seq datasets

This file contains examples of parameters and R session informations for the download and the analysis of TCGA RNA-seq datasets.

## Download datasets

This section contains an example of parameter section to use in the `TCGA-downloadDatasets.R` script to download TCGA RNA-seq datasets.

### Parameters

```
##############
# Parameters #
##############
# list of datasets
datasets <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-THCA")
# data category and type
data_category <- "Transcriptome Profiling"
data_type <- "Gene Expression Quantification"
# type of workflow
workflow_type <- "Counts"
# phenotypes of interest
phenotypes <- c("TP", "NT")
# path to output directory
output_dir <- "./output/TCGA/00-Data"
```

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
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SummarizedExperiment_1.10.1 DelayedArray_0.6.1          BiocParallel_1.14.1        
 [4] matrixStats_0.53.1          Biobase_2.40.0              GenomicRanges_1.32.3       
 [7] GenomeInfoDb_1.16.0         IRanges_2.14.10             S4Vectors_0.18.3           
[10] BiocGenerics_0.30.0         TCGAbiolinks_2.12.6        

loaded via a namespace (and not attached):
  [1] colorspace_1.3-2            ggsignif_0.6.0              selectr_0.4-1              
  [4] rjson_0.2.20                hwriter_1.3.2               circlize_0.4.4             
  [7] XVector_0.20.0              GlobalOptions_0.1.0         rstudioapi_0.7             
 [10] ggpubr_0.2.3                matlab_1.0.2                ggrepel_0.8.0              
 [13] bit64_0.9-7                 AnnotationDbi_1.46.0        fansi_0.4.0                
 [16] xml2_1.2.0                  codetools_0.2-16            splines_3.6.2              
 [19] R.methodsS3_1.7.1           doParallel_1.0.11           DESeq_1.36.0               
 [22] geneplotter_1.58.0          knitr_1.37                  jsonlite_1.7.1             
 [25] Rsamtools_1.32.2            km.ci_0.5-2                 broom_0.5.2                
 [28] annotate_1.58.0             cluster_2.1.0               R.oo_1.22.0                
 [31] readr_1.1.1                 compiler_3.6.2              httr_1.3.1                 
 [34] backports_1.1.2             assertthat_0.2.0            Matrix_1.2-17              
 [37] limma_3.36.2                cli_3.4.1                   prettyunits_1.0.2          
 [40] tools_3.6.2                 gtable_0.2.0                glue_1.6.2                 
 [43] GenomeInfoDbData_1.1.0      dplyr_1.0.10                ggthemes_3.5.0             
 [46] ShortRead_1.38.0            Rcpp_1.0.3                  vctrs_0.5.0                
 [49] Biostrings_2.48.0           nlme_3.1-142                rtracklayer_1.40.3         
 [52] iterators_1.0.9             xfun_0.29                   stringr_1.3.1              
 [55] rvest_0.3.2                 lifecycle_1.0.3             XML_3.98-1.11              
 [58] edgeR_3.22.3                zoo_1.8-6                   zlibbioc_1.26.0            
 [61] scales_0.5.0                aroma.light_3.10.0          hms_0.4.2                  
 [64] RColorBrewer_1.1-2          ComplexHeatmap_1.18.1       memoise_1.1.0              
 [67] gridExtra_2.3               KMsurv_0.1-5                ggplot2_3.3.4              
 [70] downloader_0.4              biomaRt_2.40.3              latticeExtra_0.6-28        
 [73] stringi_1.2.3               RSQLite_2.1.1               genefilter_1.66.0          
 [76] foreach_1.4.4               GenomicFeatures_1.36.4      shape_1.4.4                
 [79] rlang_1.0.6                 pkgconfig_2.0.1             bitops_1.0-6               
 [82] lattice_0.20-38             purrr_0.2.5                 GenomicAlignments_1.16.0   
 [85] bit_1.1-14                  tidyselect_1.2.0            plyr_1.8.4                 
 [88] magrittr_1.5                R6_2.2.2                    generics_0.0.2             
 [91] DBI_1.0.0                   mgcv_1.8-28                 pillar_1.8.1               
 [94] survival_3.1-7              RCurl_1.95-4.10             tibble_2.1.3               
 [97] EDASeq_2.18.0               crayon_1.3.4                survMisc_0.5.5             
[100] utf8_1.1.4                  GetoptLong_0.1.7            progress_1.2.0             
[103] locfit_1.5-9.1              grid_3.6.2                  sva_3.32.1                 
[106] data.table_1.11.4           blob_1.1.1                  ConsensusClusterPlus_1.44.0
[109] digest_0.6.15               xtable_1.8-2                tidyr_0.8.1                
[112] R.utils_2.6.0               munsell_0.5.0               survminer_0.4.6            
```

## Differential dispersion analysis

This section contains an example of parameter section to use in the `simulations-DD_analysis.R` script to identify differentially dispersed genes in a TCGA RNA-seq dataset. The parameter setting enables the anlalysis of a dataset downloaded with the parameter setting of the previous section. All datasets are available in this GitHub repository.

### Parameters

```
##############
# Parameters #
##############
# dataset
dataset <- "TCGA-BRCA"
# normalization method
normalization_method <- "TMM"
# minimum expression threshold in CPM to filter out lowly expressed genes
filter_threshold <- 1
# MDSeq parameters
## perform MDSeq outlier removal function
outlier_removal <- TRUE
## fold-change threshold to identify DD genes (not in log2 scale)
DD_FC_threshold <- 1
## p-value threshold
pval_threshold <- 0.05
## number of cores to use for computations
cores <- 4
## minimum sample size
min_sample_size <- 1
# path to output directory
output_dir <- "./output/TCGA"
```

The name of each TCGA dataset must be provided separately using the 'dataset' variable.

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
 [1] parallel  splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DiffDist_0.1.0.9000 gamlss_5.4-10       nlme_3.1-142        gamlss.dist_6.0-5   MASS_7.3-50        
 [6] gamlss.data_6.0-2   DiPhiSeq_0.2.0      MDSeq_1.0.5         lawstat_3.2         VGAM_1.0-5         
[11] mvtnorm_1.0-8       Kendall_2.2.1       edgeR_3.22.3        limma_3.36.2       

loaded via a namespace (and not attached):
[1] locfit_1.5-9.1  Rcpp_1.0.3      lattice_0.20-38 grid_3.6.2      Matrix_1.2-17   boot_1.3-28    
[7] tools_3.6.2     survival_3.1-7  compiler_3.6.2 
```
