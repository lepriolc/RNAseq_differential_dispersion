library(edgeR)
library(DiPhiSeq)
library(MDSeq)
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
## number of cores to use for computations
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
conditions <- as.factor(sim_data@sample.annotations$condition)
name <- sim_data@info.parameters$dataset
sample_size <- sim_data@info.parameters$samples.per.cond
replicate <- sim_data@info.parameters$repl.id
annotations <- sim_data@variable.annotations
## identify whether the simulated dataset is composed of highly DE genes or only lowly DE genes
highly_DE <- ifelse(sum(abs(log(annotations[, sprintf("%s.S2", mean_colname)]/annotations[, sprintf("%s.S1", mean_colname)], 2)) > log(DE_FC_threshold, 2)) > 0, TRUE, FALSE)
output_subdir_name <- ifelse(highly_DE, "00-highly_DE", "10-lowly_DE")

# filter lowly expressed genes
count_matrix <- filter_lowly_expressed_genes(sim_data@count.matrix, as.factor(conditions), normalization_method, filter_threshold)

# DiPhiSeq
diphiseq_output_dir <- sprintf("%s/10-DiPhiSeq/%s/base_effect_%s/%s/repl_%d", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (! dir.exists(diphiseq_output_dir)) {
  dir.create(diphiseq_output_dir, recursive=TRUE, mode="0775")
}
## DD gene identification
diphiseq_results <- run_DiPhiSeq(count_matrix, conditions, normalization_method, filter_threshold)
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
write.csv(diphiseq_results, file=sprintf("%s/%s_results_performance.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
diphiseq_results_perf_stats_df <- performance_stats(diphiseq_results, category_colname)
write.csv(diphiseq_results_perf_stats_df, file=sprintf("%s/%s_results_performance_stats.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=FALSE)
### AUC
diphiseq_results_auc_stats_df <- auc_stats(diphiseq_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(diphiseq_results_auc_stats_df, file=sprintf("%s/%s_results_auc_stats.csv", diphiseq_output_dir, output_basename), quote=FALSE, row.names=FALSE)

# MDSeq
## entire dataset
### normalization
exp.normalized <- normalize.counts(count_matrix, group=conditions, method=normalization_method)
mdseq_output_dir <- sprintf("%s/20-MDSeq/%s/base_effect_%s/%s/repl_%d", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
if (highly_DE) {
  mdseq_output_dir <- sprintf("%s/00-all_genes", mdseq_output_dir)
}
if (! dir.exists(mdseq_output_dir)) {
  dir.create(mdseq_output_dir, recursive=TRUE, mode="0775")
}
### fit MDSeq model
output_basename <- sprintf("%s_%s_MDSeq", name, normalization_method)
fit <- MDSeq_fit(exp.normalized, conditions, outlier_removal, NULL, NULL, cores, min_sample_size, mdseq_output_dir, output_basename)
write.csv(fit$Dat, file=sprintf("%s/%s_fit_Dat.csv", mdseq_output_dir, output_basename), quote=FALSE, row.names=TRUE)
### DD gene identification
mdseq_results <- MDSeq_DE_DD(fit, conditions, DD_FC_threshold)
write.csv(mdseq_results, file=sprintf("%s/%s_FC_%s.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)
### performance
DD_pval_colanme <- "FDR.dispersion" # use p-values corrected by the Benjamini-Yekutieli FDR-controlling procedure
#### DD analysis results with contingency categories
mdseq_results <- add_DD_performance_to_results(mdseq_results, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
write.csv(mdseq_results, file=sprintf("%s/%s_FC_%s_performance.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)
#### DD analysis performance statistics
category_colname <- sprintf("category.%s", DD_pval_colanme)
mdseq_results_perf_stats_df <- performance_stats(mdseq_results, category_colname)
write.csv(mdseq_results_perf_stats_df, file=sprintf("%s/%s_FC_%s_performance_stats.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=FALSE)
#### AUC
mdseq_results_auc_stats_df <- auc_stats(mdseq_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
write.csv(mdseq_results_auc_stats_df, file=sprintf("%s/%s_FC_%s_auc_stats.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=FALSE)

## DD analysis for lowly DE genes
if (highly_DE) {
  mdseq_output_dir <- sprintf("%s/20-MDSeq/%s/base_effect_%s/%s/repl_%d/10-lowly_DE_genes", output_dir, output_subdir_name, sub("[.]", "_", DE_FC_threshold), sample_size, replicate)
  output_basename <- sprintf("%s_lowly_DE", output_basename)
  if (! dir.exists(mdseq_output_dir)) {
    dir.create(mdseq_output_dir, recursive=TRUE, mode="0775")
  }
  ### normalization, fit MDSeq model and DD gene identification
  mdseq_DD_lowly_DE_results <-MDSeq_DD_for_lowly_DE(count_matrix, conditions, normalization_method, outlier_removal, NULL, NULL, cores, min_sample_size, DE_FC_threshold, DD_FC_threshold, pval_threshold, mdseq_output_dir, name)
  write.csv(mdseq_DD_lowly_DE_results, file=sprintf("%s/%s_FC_%s.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)
  ### performance
  DD_pval_colanme <- "FDR.dispersion"
  #### DD analysis results with contingency categories
  mdseq_DD_lowly_DE_results <- add_DD_performance_to_results(mdseq_DD_lowly_DE_results, annotations, DD_FC_threshold, pval_threshold, dispersion_colname, DD_pval_colanme)
  write.csv(mdseq_DD_lowly_DE_results, file=sprintf("%s/%s_FC_%s_performance.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=TRUE)
  #### DD analysis performance statistics
  category_colname <- sprintf("category.%s", DD_pval_colanme)
  mdseq_DD_lowly_DE_results_perf_stats_df <- performance_stats(mdseq_DD_lowly_DE_results, category_colname)
  write.csv(mdseq_DD_lowly_DE_results_perf_stats_df, file=sprintf("%s/%s_FC_%s_performance_stats.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=FALSE)
  #### AUC
  mdseq_DD_lowly_DE_results_auc_stats_df <- auc_stats(mdseq_DD_lowly_DE_results, DD_pval_colanme, "labels.DD", annotations, mean_colname, DE_FC_threshold)
  write.csv(mdseq_DD_lowly_DE_results_auc_stats_df, file=sprintf("%s/%s_FC_%s_auc_stats.csv", mdseq_output_dir, output_basename, sub("[.]", "_", DD_FC_threshold)), quote=FALSE, row.names=FALSE)
}

