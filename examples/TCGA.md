# Examples of parameters for the analysis of TCGA RNA-seq datasets

This file contains examples of parameters and R session informations for the download and the analysis of TCGA RNA-seq datasets.

## Download datasets

This section contains an example of parameter section to use in the analysis script `TCGA-downloadDatasets.R` and the R session information to download TCGA RNA-seq datasets.

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
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8     LC_MONETARY=fr_FR.UTF-8   
 [6] LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SummarizedExperiment_1.10.1 DelayedArray_0.6.1          BiocParallel_1.14.1         matrixStats_0.53.1         
 [5] Biobase_2.40.0              GenomicRanges_1.32.3        GenomeInfoDb_1.16.0         IRanges_2.14.10            
 [9] S4Vectors_0.18.3            BiocGenerics_0.30.0         TCGAbiolinks_2.12.6        

loaded via a namespace (and not attached):
  [1] colorspace_1.3-2            selectr_0.4-1               ggsignif_0.6.0              rjson_0.2.20               
  [5] hwriter_1.3.2               circlize_0.4.4              XVector_0.20.0              GlobalOptions_0.1.0        
  [9] rstudioapi_0.7              ggpubr_0.2.3                matlab_1.0.2                ggrepel_0.8.0              
 [13] bit64_0.9-7                 AnnotationDbi_1.46.0        xml2_1.2.0                  codetools_0.2-16           
 [17] splines_3.6.2               R.methodsS3_1.7.1           doParallel_1.0.11           DESeq_1.36.0               
 [21] geneplotter_1.58.0          knitr_1.37                  jsonlite_1.7.1              Rsamtools_1.32.2           
 [25] km.ci_0.5-2                 broom_0.5.2                 annotate_1.58.0             cluster_2.1.0              
 [29] R.oo_1.22.0                 readr_1.1.1                 compiler_3.6.2              httr_1.3.1                 
 [33] backports_1.1.2             assertthat_0.2.0            Matrix_1.2-17               limma_3.36.2               
 [37] prettyunits_1.0.2           tools_3.6.2                 gtable_0.2.0                glue_1.3.1                 
 [41] GenomeInfoDbData_1.1.0      dplyr_0.8.3                 ggthemes_3.5.0              ShortRead_1.38.0           
 [45] Rcpp_1.0.3                  Biostrings_2.48.0           nlme_3.1-142                rtracklayer_1.40.3         
 [49] iterators_1.0.9             xfun_0.29                   stringr_1.3.1               rvest_0.3.2                
 [53] XML_3.98-1.11               edgeR_3.22.3                zoo_1.8-6                   zlibbioc_1.26.0            
 [57] scales_0.5.0                aroma.light_3.10.0          hms_0.4.2                   RColorBrewer_1.1-2         
 [61] ComplexHeatmap_1.18.1       memoise_1.1.0               gridExtra_2.3               KMsurv_0.1-5               
 [65] ggplot2_3.3.4               downloader_0.4              biomaRt_2.40.3              latticeExtra_0.6-28        
 [69] stringi_1.2.3               RSQLite_2.1.1               genefilter_1.66.0           foreach_1.4.4              
 [73] GenomicFeatures_1.36.4      shape_1.4.4                 rlang_0.4.11                pkgconfig_2.0.1            
 [77] bitops_1.0-6                lattice_0.20-38             purrr_0.2.5                 GenomicAlignments_1.16.0   
 [81] bit_1.1-14                  tidyselect_0.2.5            plyr_1.8.4                  magrittr_1.5               
 [85] R6_2.2.2                    generics_0.0.2              DBI_1.0.0                   mgcv_1.8-28                
 [89] pillar_1.4.2                survival_3.1-7              RCurl_1.95-4.10             tibble_2.1.3               
 [93] EDASeq_2.18.0               crayon_1.3.4                survMisc_0.5.5              GetoptLong_0.1.7           
 [97] progress_1.2.0              locfit_1.5-9.1              grid_3.6.2                  sva_3.32.1                 
[101] data.table_1.11.4           blob_1.1.1                  ConsensusClusterPlus_1.44.0 digest_0.6.15              
[105] xtable_1.8-2                tidyr_0.8.1                 R.utils_2.6.0               munsell_0.5.0              
[109] survminer_0.4.6
```

## Differential dispersion analysis

This section contains an example of parameter section to use in the analysis script `simulations-DD_analysis.R` and the R session information to identify differentially dispersed genes in TCGA-BRCA RNA-seq dataset.

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
## fold-change threshold to identify highly genes (not in log2 scale)
DE_FC_threshold <- ifelse(dataset=="TCGA-BRCA", 1.25, 1.3)
## number of cores to use for computations
cores <- 4
## minimum sample size
min_sample_size <- 1
# path to output directory
output_dir <- "./output/TCGA"
```

The name of each TCGA dataset must be provided separately using the 'dataset' variable. The datasets analyzed in the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression' are: TCGA-BRCA, TCGA-COAD, TCGA-HNSC, TCGA-KIRC, TCGA-KIRP, TCGA-LIHC, TCGA-LUAD, TCGA-LUSC, TCGA-PRAD, TCGA-THCA.

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
[1] MDSeq_1.0.5    DiPhiSeq_0.2.0 edgeR_3.22.3   limma_3.36.2  

loaded via a namespace (and not attached):
[1] compiler_3.6.2  tools_3.6.2     Rcpp_1.0.3      grid_3.6.2      locfit_1.5-9.1  lattice_0.20-38
```

