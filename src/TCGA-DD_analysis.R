library(edgeR)
library(DiPhiSeq)
library(MDSeq)


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


#############
# Functions #
#############
source("./src/DD_analysis-functions.R")


#############
# Load data #
#############

# expression data
dataset_input_dir <- sprintf("%s/00-Data/%s", output_dir, dataset)
exp_Rdata_file <- sprintf("%s/%s_Counts_exp_pheno-patients_TP_NT_samples.RData", dataset_input_dir, dataset)
load(exp_Rdata_file) # contains 'exp_matrix' and 'metadata_pheno_df'

# output basename
pheno_1 <- levels(metadata_pheno_df$phenotype)[1]
pheno_2 <- levels(metadata_pheno_df$phenotype)[2]
output_basename <- sprintf("%s_%s_vs_%s", dataset, pheno_1, pheno_2)


############
# Analysis #
############
# remove genes with read count = 0 for all samples
exp_matrix_no0row <- exp_matrix[!apply(exp_matrix, 1, function(x) { sum(x) == 0 }),]

# filter lowly expressed genes
conditions <- metadata_pheno_df$phenotype
count_matrix <- filter_lowly_expressed_genes(exp_matrix_no0row, conditions, normalization_method, filter_threshold)


# DiPhiSeq
## output directory
diphiseq_output_dir <- sprintf("%s/10-DiPhiSeq/%s", output_dir, dataset)
if (! dir.exists(diphiseq_output_dir)) {
  dir.create(diphiseq_output_dir, recursive=TRUE, mode="0775")
}
## conditions for DiPhiSeq: a numeric vector of 1s and 2s
diphiseq_conditions <- as.character(conditions)
pheno_1 <- levels(conditions)[1]
pheno_2 <- levels(conditions)[2]
diphiseq_conditions[which(diphiseq_conditions==pheno_1)] <- 1
diphiseq_conditions[which(diphiseq_conditions==pheno_2)] <- 2
diphiseq_conditions <- as.numeric(diphiseq_conditions)
## DD gene identification
diphiseq_results <- run_DiPhiSeq(count_matrix, diphiseq_conditions, normalization_method, filter_threshold)
### rename columns in tab
colnames(diphiseq_results$results$tab)[which(colnames(diphiseq_results$results$tab)=="phi1")] <- sprintf("phi.%s", pheno_1)
colnames(diphiseq_results$results$tab)[which(colnames(diphiseq_results$results$tab)=="phi2")] <- sprintf("phi.%s", pheno_2)
colnames(diphiseq_results$results$tab)[which(colnames(diphiseq_results$results$tab)=="beta1")] <- sprintf("beta.%s", pheno_1)
colnames(diphiseq_results$results$tab)[which(colnames(diphiseq_results$results$tab)=="beta2")] <- sprintf("beta.%s", pheno_2)
## write outputs
if (identical(rownames(diphiseq_results$results$mumat), rownames(diphiseq_results$results$phimat))) {
  if (dim(diphiseq_results$results$tab)[1] == dim(diphiseq_results$results$mumat)[1]) {
    rownames(diphiseq_results$results$tab) <- rownames(diphiseq_results$results$mumat)
  } else {
    stop("DiPhiSeq tab and mu matrices have different number of rows !")
  }
} else {
  stop("Rownames of DiPhiSeq mu and phi matrices are different !")
}
output_basename <- sprintf("%s_%s_vs_%s_%s_DiPhiSeq", dataset, pheno_1, pheno_2, normalization_method)
write.csv(diphiseq_results$results$tab, file=sprintf("%s/%s_results.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
write.csv(diphiseq_results$results$mumat, file=sprintf("%s/%s_mumat.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
write.csv(diphiseq_results$results$phimat, file=sprintf("%s/%s_phimat.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
write.csv(diphiseq_results$samples, file=sprintf("%s/%s_samples.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)


# MDSeq
## output directory
mdseq_output_dir <- sprintf("%s/20-MDSeq/%s", output_dir, dataset)
if (! dir.exists(mdseq_output_dir)) {
  dir.create(mdseq_output_dir, recursive=TRUE, mode="0775")
}
## write metadata
write.csv(metadata_pheno_df, file=sprintf("%s/%s_metadata.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)

## design matrix
design <- get.model.matrix(conditions)
## covariate matrix
covariate <- "batch"
### covariate in a matrix with numeric values
cov_matrix <- as.matrix(sapply(metadata_pheno_df[, covariate], as.numeric))
colnames(cov_matrix) <- covariate
rownames(cov_matrix) <- rownames(metadata_pheno_df)
write.csv(cov_matrix, file=sprintf("%s/%s_MDSeq_cov_matrix_numeric.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)

## DD analysis
### entire datasets
analysis_output_dir <- sprintf("%s/00-all_genes", mdseq_output_dir)
if (! dir.exists(analysis_output_dir)) {
  dir.create(analysis_output_dir, recursive=TRUE, mode="0775")
}
analysis_output_basename <- sprintf("%s_%s_MDSeq", output_basename, normalization_method)
#### normalization
exp.normalized <- normalize.counts(count_matrix, group=conditions, method=normalization_method)
#### fit MDSeq model
fit <- MDSeq_fit(exp.normalized, conditions, outlier_removal, cov_matrix, cov_matrix, cores, min_sample_size, analysis_output_dir, analysis_output_basename)
write.csv(fit$Dat, file=sprintf("%s/%s_fit_Dat.csv", analysis_output_dir, analysis_output_basename), quote=FALSE, row.names=TRUE)
#### DD gene identification
mdseq_results <- MDSeq_DE_DD(fit, conditions, DD_FC_threshold)
write.csv(mdseq_results, file=sprintf("%s/%s_FC_%s.csv", analysis_output_dir, analysis_output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)

### lowly DE genes
analysis_output_dir <- sprintf("%s/10-lowly_DE_genes", mdseq_output_dir)
if (! dir.exists(analysis_output_dir)) {
  dir.create(analysis_output_dir, recursive=TRUE, mode="0775")
}
#### normalization, fit MDSeq model and DD gene identification
mdseq_DD_lowly_DE_results <-MDSeq_DD_for_lowly_DE(count_matrix, conditions, normalization_method, outlier_removal, cov_matrix, cov_matrix, cores, min_sample_size, DE_FC_threshold, DD_FC_threshold, pval_threshold, analysis_output_dir, output_basename)
write.csv(mdseq_DD_lowly_DE_results, file=sprintf("%s/%s_%s_MDSeq_lowly_DE_FC_%s.csv", analysis_output_dir, output_basename, normalization_method, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)

