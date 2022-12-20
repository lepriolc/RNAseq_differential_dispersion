library(edgeR)
library(lawstat)
library(MDSeq)
library(DiPhiSeq)
library(gamlss)
library(DiffDist)


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

# conditions
conditions <- metadata_pheno_df$phenotype
pheno_1 <- levels(conditions)[1]
pheno_2 <- levels(conditions)[2]


############
# Analysis #
############
# remove genes with read count = 0 for all samples
exp_matrix_no0row <- exp_matrix[!apply(exp_matrix, 1, function(x) { sum(x) == 0 }),]

# filter lowly expressed genes
count_matrix <- filter_lowly_expressed_genes(exp_matrix_no0row, conditions, normalization_method, filter_threshold)


# Levene's test
## output directory
Levene_output_dir <- sprintf("%s/10-Levene/%s", output_dir, dataset)
if (! dir.exists(Levene_output_dir)) {
  dir.create(Levene_output_dir, recursive=TRUE, mode="0775")
}
## DD gene identification
condition_1_samples <- rownames(metadata_pheno_df[which(metadata_pheno_df$phenotype==pheno_1),])
condition_2_samples <- rownames(metadata_pheno_df[which(metadata_pheno_df$phenotype==pheno_2),])
Levene_results <- run_Levene(count_matrix, conditions, condition_1_samples, condition_2_samples, normalization_method)
## write outputs
output_basename <- sprintf("%s_%s_vs_%s_%s_Levene_test", dataset, pheno_1, pheno_2, normalization_method)
write.csv(Levene_results, file=sprintf("%s/%s_results.csv", Levene_output_dir, output_basename), quote=FALSE, row.names=TRUE)


# MDSeq
## output directory
mdseq_output_dir <- sprintf("%s/20-MDSeq/%s", output_dir, dataset)
if (! dir.exists(mdseq_output_dir)) {
  dir.create(mdseq_output_dir, recursive=TRUE, mode="0775")
}
## write metadata
output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq", dataset, pheno_1, pheno_2, normalization_method)
write.csv(metadata_pheno_df, file=sprintf("%s/%s_metadata.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)

## design matrix
design <- get.model.matrix(conditions)
## covariate matrix
covariate <- "batch"
### covariate in a matrix with numeric values
cov_matrix <- as.matrix(sapply(metadata_pheno_df[, covariate], as.numeric))
colnames(cov_matrix) <- covariate
rownames(cov_matrix) <- rownames(metadata_pheno_df)
write.csv(cov_matrix, file=sprintf("%s/%s_cov_matrix_numeric.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)

## DD analysis
### normalization
exp.normalized <- normalize.counts(count_matrix, group=conditions, method=normalization_method)
### fit MDSeq model
fit <- run_MDSeq_fit(exp.normalized, conditions, outlier_removal, cov_matrix, cov_matrix, cores, min_sample_size, mdseq_output_dir, output_basename)
write.csv(fit$Dat, file=sprintf("%s/%s_fit.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD gene identification
mdseq_results <- run_MDSeq_DD(fit, conditions, DD_FC_threshold)
write.csv(mdseq_results, file=sprintf("%s/%s_results.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)


# DiPhiSeq
count_matrix <- count_matrix[1:200,]
## output directory
diphiseq_output_dir <- sprintf("%s/30-DiPhiSeq/%s", output_dir, dataset)
if (! dir.exists(diphiseq_output_dir)) {
  dir.create(diphiseq_output_dir, recursive=TRUE, mode="0775")
}
## conditions for DiPhiSeq: a numeric vector of 1s and 2s
diphiseq_conditions <- as.character(conditions)
diphiseq_conditions[which(diphiseq_conditions==pheno_1)] <- 1
diphiseq_conditions[which(diphiseq_conditions==pheno_2)] <- 2
diphiseq_conditions <- as.numeric(diphiseq_conditions)
## DD gene identification
diphiseq_results <- run_DiPhiSeq(count_matrix, diphiseq_conditions, normalization_method, filter_threshold) # unsused filter_threshold parameter
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


# GAMLSS
## output directory
gamlss_output_dir <- sprintf("%s/40-GAMLSS/%s", output_dir, dataset)
if (! dir.exists(gamlss_output_dir)) {
  dir.create(gamlss_output_dir, recursive=TRUE, mode="0775")
}
## DD gene identification
gamlss_results <- run_GAMLSS(count_matrix, conditions, pheno_2, normalization_method)
## write outputs
output_basename <- sprintf("%s_%s_vs_%s_%s_GAMLSS", dataset, pheno_1, pheno_2, normalization_method)
write.csv(gamlss_results, file=sprintf("%s/%s_results.csv", gamlss_output_dir, output_basename), quote=FALSE, row.names=TRUE)


# DiffDist
## output directory
diffdist_output_dir <- sprintf("%s/50-DiffDist/%s", output_dir, dataset)
if (! dir.exists(diffdist_output_dir)) {
  dir.create(diffdist_output_dir, recursive=TRUE, mode="0775")
}
## conditions for DiffDist: a numeric vector of 1s and 2s
diffdist_conditions <- as.character(conditions)
diffdist_conditions[which(diffdist_conditions==pheno_1)] <- 1
diffdist_conditions[which(diffdist_conditions==pheno_2)] <- 2
diffdist_conditions <- as.factor(diffdist_conditions)
## DD gene identification
diffdist_results <- run_DiffDist(count_matrix, diffdist_conditions, normalization_method)
## replace numeric conditions by conditions in column names
colnames(diffdist_results)[which(grepl("1vs2", colnames(diffdist_results)))] <- sub("1vs2", sprintf("%svs%s", pheno_1, pheno_2), colnames(diffdist_results)[which(grepl("1vs2", colnames(diffdist_results)))])
## write outputs
output_basename <- sprintf("%s_%s_vs_%s_%s_DiffDist", dataset, pheno_1, pheno_2, normalization_method)
write.csv(diffdist_results, file=sprintf("%s/%s_results.csv", diffdist_output_dir, output_basename), quote=FALSE, row.names=TRUE)

