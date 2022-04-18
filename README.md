# RNA-seq differential dispersion

## Content

This repository contains all the scripts to reproduce the results published in the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression'.

## Code

All the R scripts are stored in the src/ directory.

There are 3 types of scripts:

- function files
- analysis scripts
- figure scripts

### Function files

The function files only contain the defintion of functions called by the analysis and figures scripts. Their names end by '-functions.R'.

The file `DD_analysis-functions.R` contains generic functions for the identification of differentially dispersed (DD) genes in RNA-seq datasets. The analysis of simulated and TCGA datasets requires this file. The `simulations-functions.R` file is dedicated to the analysis of simulated datasets.

### Analysis and figures scripts

The analysis and figures R scripts contain different sections. The 'Parameter' section defines the parameters defined for the analysis. The 'Analysis' section contains the commands to run the analysis or generate the figures published in the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression' and is not supposed to be modified to reproduce them.

The names of anlalysis scripts end by '-analysis.R' and the names of figure scripts end by '-figures.R'

#### Simulated datasets

By default, all the files generated by the analysis and figures scripts are stored in the directory './output/simulations/'.

The script `simulations-generateDatasets.R` generates simulated RNA-seq datasets based on the parameters provided in its 'Parameters' section. The outputs are stored in the '00-Data/' subdirectory.

The script `simulations-DD_analysis.R` identifies DD genes in simulated RNA-seq datasets using DiPhiSeq and MDSeq and evaluates the performances of these methods. It requires the functions defined in the DD_analysis-functions.R and simulations-functions.R files and parameters provided in the 'Parameters' section. The outputs of DiPhiSeq and MDSeq are stored in the subdirectories '10-DiPhiSeq/' and '20-MDSeq/' respectively, one subdirectory per dataset. The subdirectory structure is set according to the type of datasets (either containing highly differentially expressed (DE) genes or only lowly DE genes) and the parameters used.

#### TCGA datasets

By default, all the files generated by the analysis and figures scripts are stored in the directory './output/TCGA/'.

The script `TCGA-downloadDatasets.R` downloads TCGA RNA-seq datasets based on the parameters provided in its 'Parameters' section. The downloaded files are stored in RData files in the '00-Data/' subdirectory. These files, which are the inputs of the differential dispersion analysis, are available in this GitHub repository.

The file `TCGA-DD_analysis.R` identifies DD genes in a TCGA RNA-seq dataset using DiPhiSeq and MDSeq. The name of the TCGA dataset is defined by the 'dataset' variable and the input files are the outputs files of the TCGA-downloadDatasets.R script stored in the directory whose path is defined in the variable 'dataset_input_dir'. The outputs of DiPhiSeq and MDSeq are stored in the subdirectories '10-DiPhiSeq/' and '20-MDSeq/' respectively, one subdirectory per dataset.

The file `TCGA-GO_cluster.R` performs Gene Ontology term enrichment analysis for overdispersed (DD+) and underdispersed genes in tumors (DD-) among lowly differentially expressed (DE) genes, highly upregulated (DE+) and highly downregulated genes in tumors (DE-) respectively for all TCGA datasets defined in the 'datasets' variable. To ease comparison, representative terms of closely related GO terms are identified thanks to semantic similarity and hierarchical clustering. This script requires parameters provided its 'Parameters' section and functions defined in the TCGA-functions.R file. The input files are the DiPhiSeq and MDSeq output files generated by the TCGA-DD_analysis.R script and stored in subdirectories of the directory whose path is defined in the variable 'output_dir'. The outputs are stored in the subdirectory '30-GO_cluster/'. The subdirectory '00-gene-lists/' contains the identifiers of genes identified either by DiPhiSeq or MDSeq in categories according differential expression or dispersion. The subdirectory '10-GO_analysis/' contains the outputs of GO term cluster analysis for DD+, DD-, DE+ and DE- genes. The first page of '30-GO_cluster/10-GO_analysis/lowly_DE_DD+/20-GO_cluster/TCGA_TP_vs_NT_FC_1_enrichGO_impified_Rel_0_8_cluster_generic_terms_customplot.pdf' file is the figure 6 of the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression'.

The script `TCGA-figures.R` generates the figures 4, 5, S4, S5, S6 and S7 based on the the DiPhiSeq and MDSeq output files generated by the TCGA-DD_analysis.R script. The figures are contained in files stored in the directory whose path is defined in the variable 'output_dir'.

## Examples

The examples/ directory contains the parameter sections of scripts for the analysis of a simulated dataset and a TCGA dataset.

## Contact

For more details about the analysis, please refer to the article entitled 'Detection of genes with a differential expression dispersion unravels the role of autophagy in cancer progression'

If you have any question regarding the scripts in this repository, please contact chris.lepriol@gmail.com.

