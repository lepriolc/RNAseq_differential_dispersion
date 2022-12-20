library(compcodeR)
library(edgeR)
library(lawstat)
library(MDSeq)
library(DiPhiSeq)
library(gamlss)
library(DiffDist)
library(ROCR)


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
## number of cores to use for MDSeq computations
cores <- 4
## minimum sample size
min_sample_size <- 1
# path to output directory
output_dir <- "./output/simulations"


#############
# Functions #
#############
source("./src/DD_analysis-functions.R")
source("./src/simulations-functions.R")


############
# Analysis #
############

# load simulated dataset and get simulation parameters
sim_data <- readRDS(rds_file) # contains compData object (compcodeR)
sample_annotations <- sim_data@sample.annotations
conditions <- as.factor(sample_annotations$condition)
name <- sim_data@info.parameters$dataset
sample_size <- sim_data@info.parameters$samples.per.cond
replicate <- sim_data@info.parameters$repl.id
annotations <- sim_data@variable.annotations
## identify whether the simulated dataset is composed of highly DE genes or only lowly DE genes
highly_DE <- ifelse(sum(abs(log(annotations[, sprintf("%s.S2", mean_colname)]/annotations[, sprintf("%s.S1", mean_colname)], 2)) > log(DE_FC_threshold, 2)) > 0, TRUE, FALSE)
output_subdir_name <- ifelse(highly_DE, "00-Highly_DE", "10-Lowly_DE")

# filter lowly expressed genes
count_matrix <- filter_lowly_expressed_genes(sim_data@count.matrix, as.factor(conditions), normalization_method, filter_threshold)


# Levene's test
Levene_output_dir <- sprintf("%s/10-Levene/%s/base_effect_%s/%s/repl_%d", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (! dir.exists(Levene_output_dir)) {
  dir.create(Levene_output_dir, recursive=TRUE, mode="0775")
}
## run Levene's test
condition_1_samples <- rownames(sample_annotations[which(sample_annotations$condition==levels(conditions)[1]),])
condition_2_samples <- rownames(sample_annotations[which(sample_annotations$condition==levels(conditions)[2]),])
Levene_results <- run_Levene(count_matrix, conditions, condition_1_samples, condition_2_samples, normalization_method)
## write outputs
output_basename <- sprintf("%s_%s_Levene_test_results", name, normalization_method)
write.csv(Levene_results, file=sprintf("%s/%s.csv", Levene_output_dir, output_basename), quote=FALSE, row.names=TRUE)
## performance
DD_pval_colanme <- "p_value.BH" # use p-values corrected by the Benjamini-Hochberg FDR-controlling procedure
### DD analysis results with performance
Levene_results <- add_DD_performance_to_results(Levene_results, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
write.csv(Levene_results, file=sprintf("%s/%s_performance.csv", Levene_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
Levene_results_perf_stats_df <- performance_stats(Levene_results, category_colname)
write.csv(Levene_results_perf_stats_df, file=sprintf("%s/%s_performance_stats.csv", Levene_output_dir, output_basename), quote=FALSE, row.names=FALSE)
### AUC
Levene_results_auc_stats_df <- auc_stats(Levene_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(Levene_results_auc_stats_df, file=sprintf("%s/%s_auc_stats.csv", Levene_output_dir, output_basename), quote=FALSE, row.names=FALSE)


# MDSeq
## normalization
exp.normalized <- normalize.counts(count_matrix, group=conditions, method=normalization_method)
### warning: normalize counts with MDSeq normalize.counts may be slightly different from those obtained with edgeR approach
mdseq_output_dir <- sprintf("%s/20-MDSeq/%s/base_effect_%s/%s/repl_%d", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (! dir.exists(mdseq_output_dir)) {
  dir.create(mdseq_output_dir, recursive=TRUE, mode="0775")
}
## fit MDSeq model
output_basename <- sprintf("%s_%s_MDSeq", name, normalization_method)
fit <- run_MDSeq_fit(exp.normalized, conditions, outlier_removal, NULL, NULL, cores, min_sample_size, mdseq_output_dir, output_basename)
write.csv(fit$Dat, file=sprintf("%s/%s_fit.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
## DD analysis with the entire dataset
output_basename <- sprintf("%s_results", output_basename)
mdseq_results <- run_MDSeq_DD(fit, conditions, DD_FC_threshold)
write.csv(mdseq_results, file=sprintf("%s/%s.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
## performance
DD_pval_colanme <- "FDR.dispersion" # use p-values corrected by the Benjamini-Yekutieli FDR-controlling procedure
### DD analysis results with contingency categories
mdseq_results <- add_DD_performance_to_results(mdseq_results, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
write.csv(mdseq_results, file=sprintf("%s/%s_performance.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
mdseq_results_perf_stats_df <- performance_stats(mdseq_results, category_colname)
write.csv(mdseq_results_perf_stats_df, file=sprintf("%s/%s_performance_stats.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=FALSE)
### AUC
mdseq_results_auc_stats_df <- auc_stats(mdseq_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(mdseq_results_auc_stats_df, file=sprintf("%s/%s_auc_stats.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=FALSE)


# DiPhiSeq
diphiseq_output_dir <- sprintf("%s/30-DiPhiSeq/%s/base_effect_%s/%s/repl_%d", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (! dir.exists(diphiseq_output_dir)) {
  dir.create(diphiseq_output_dir, recursive=TRUE, mode="0775")
}
## run DiPhiSeq
diphiseq_results <- run_DiPhiSeq(count_matrix, conditions, normalization_method, filter_threshold) # unused filter_threshold parameter
## write outputs
output_basename <- sprintf("%s_%s_DiPhiSeq", name, normalization_method)
if (identical(rownames(diphiseq_results$results$mumat), rownames(diphiseq_results$results$phimat))) {
  if (dim(diphiseq_results$results$tab)[1] == dim(diphiseq_results$results$mumat)[1]) {
    rownames(diphiseq_results$results$tab) <- rownames(diphiseq_results$results$mumat)
  } else {
    stop("DiPhiSeq tab and mu matrices have different number of rows !")
  }
} else {
  stop("Rownames of DiPhiSeq mu and phi matrices are different !")
}
write.csv(diphiseq_results$results$tab, file=sprintf("%s/%s_results.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
write.csv(diphiseq_results$results$mumat, file=sprintf("%s/%s_mumat.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
write.csv(diphiseq_results$results$phimat, file=sprintf("%s/%s_phimat.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
write.csv(diphiseq_results$samples, file=sprintf("%s/%s_samples.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
## performance
DD_pval_colanme <- "fdr.phi" # use p-values corrected by the Benjamini-Hochberg FDR-controlling procedure
### DD analysis results with performance
diphiseq_results <- add_DD_performance_to_results(diphiseq_results$results$tab, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
output_basename <- sprintf("%s_results", output_basename)
write.csv(diphiseq_results, file=sprintf("%s/%s_performance.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
diphiseq_results_perf_stats_df <- performance_stats(diphiseq_results, category_colname)
write.csv(diphiseq_results_perf_stats_df, file=sprintf("%s/%s_performance_stats.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=FALSE)
### AUC
diphiseq_results_auc_stats_df <- auc_stats(diphiseq_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(diphiseq_results_auc_stats_df, file=sprintf("%s/%s_auc_stats.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=FALSE)


# GAMLSS
gamlss_output_dir <- sprintf("%s/40-GAMLSS/%s/base_effect_%s/%s/repl_%d", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (! dir.exists(gamlss_output_dir)) {
  dir.create(gamlss_output_dir, recursive=TRUE, mode="0775")
}
## run GAMLSS
gamlss_results <- run_GAMLSS(count_matrix, conditions, levels(conditions)[1], normalization_method)
## write outputs
output_basename <- sprintf("%s_%s_GAMLSS_results", name, normalization_method)
write.csv(gamlss_results, file=sprintf("%s/%s.csv", gamlss_output_dir, output_basename), quote=FALSE, row.names=TRUE)
## performance
DD_pval_colanme <- "padj.cv" # use p-values corrected by the Benjamini-Hochberg FDR-controlling procedure
### DD analysis results with performance
gamlss_results <- add_DD_performance_to_results(gamlss_results, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
write.csv(gamlss_results, file=sprintf("%s/%s_performance.csv", gamlss_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
gamlss_results_perf_stats_df <- performance_stats(gamlss_results, category_colname)
write.csv(gamlss_results_perf_stats_df, file=sprintf("%s/%s_performance_stats.csv", gamlss_output_dir, output_basename), quote=FALSE, row.names=FALSE)
### AUC
gamlss_results_auc_stats_df <- auc_stats(gamlss_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(gamlss_results_auc_stats_df, file=sprintf("%s/%s_auc_stats.csv", gamlss_output_dir, output_basename), quote=FALSE, row.names=FALSE)


# DiffDist
diffdist_output_dir <- sprintf("%s/50-DiffDist/%s/base_effect_%s/%s/repl_%d_3", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (! dir.exists(diffdist_output_dir)) {
  dir.create(diffdist_output_dir, recursive=TRUE, mode="0775")
}
## run DiffDist
diffdist_results <- run_DiffDist(count_matrix, conditions, normalization_method)
## write outputs
output_basename <- sprintf("%s_%s_DiffDist_results", name, normalization_method)
write.csv(diffdist_results, file=sprintf("%s/%s.csv", diffdist_output_dir, output_basename), quote=FALSE, row.names=TRUE)
## performance
DD_pval_colanme <- sprintf("disp.%svs%s.pval.BH", levels(conditions)[1], levels(conditions)[2]) # use p-values corrected by the Benjamini-Hochberg FDR-controlling procedure
### DD analysis results with performance
diffdist_results <- add_DD_performance_to_results(diffdist_results, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
write.csv(diffdist_results, file=sprintf("%s/%s_performance.csv", diffdist_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
diffdist_results_perf_stats_df <- performance_stats(diffdist_results, category_colname)
write.csv(diffdist_results_perf_stats_df, file=sprintf("%s/%s_performance_stats.csv", diffdist_output_dir, output_basename), quote=FALSE, row.names=FALSE)
### AUC
diffdist_results_auc_stats_df <- auc_stats(diffdist_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(diffdist_results_auc_stats_df, file=sprintf("%s/%s_auc_stats.csv", diffdist_output_dir, output_basename), quote=FALSE, row.names=FALSE)

