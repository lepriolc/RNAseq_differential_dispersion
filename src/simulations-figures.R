library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(UpSetR)


##############
# Parameters #
##############
# simulations
## sample sizes
sample_sizes <- c(20, 30, 40, 50, 100)
## number of genes
genes <- 10000
## sequencing depth
depth <- 5000000
## DE genes
non_DE_runif_FC <- TRUE # generate FC for non DE genes using runif() function
## DD genes
DD_fraction <- 0.5
# DD_down_fraction <- 0.5 # fraction of DD with a decrease of dispersion in the second condition
DD_base_effect <- 1.5
DD_effect_param <- 1
non_DD_runif_FC <- FALSE # generate FC for non DD genes using runif() function
## outliers
single_outlier_high_fraction <- 0.1
single_outlier_low_fraction <- 0
random_outlier_high_fraction <- 0
random_outlier_low_fraction <- 0
## number of replicates
replicates <- 10

# DD analysis
## normalization method
normalization_method <- "TMM"
## minimum expression threshold in CPM to filter out lowly expressed genes
filter_threshold <- 1
## p-value threshold
pval_threshold <- 0.05
## DiPhiSeq parameters
### p-value correction
diphiseq_pval_correction <- "BH"
## MDSeq parameters
### perform MDSeq outlier removal function
outlier_removal <- TRUE
### fold-change threshold to identify DD genes (not in log2 scale)
DD_FC_threshold <- 1
### p-value correction
mdseq_pval_correction <- "BY"

# method order to display in plots
method_order <- c("Levene", "MDSeq", "DiPhiSeq", "GAMLSS", "DiffDist")

# path to output directory
output_dir <- "./output/simulations"


#############
# Functions #
#############
source("./src/simulations-functions.R")


############
# Analysis #
############
non_DE_runif_FC_name <- ifelse(non_DE_runif_FC, "runif", "1")
non_DD_runif_FC_name <- ifelse(non_DD_runif_FC, "runif", "1")
outlier_removal_bool_param <- ifelse(outlier_removal, 1, 0)

# Figure 1: Ability to identify differentially dispersed genes.
figure_basename <- "fig1"
## simulation parameters
### DE genes
DE_fraction <- 0.5

pdf(sprintf("%s/%s.pdf", output_dir, figure_basename), width=9)
for (DE_base_effect in c(1.1, 1.2, 1.3, 1.4, 1.5)) {
  DE_effect_param <- get_DE_effect_param(DE_base_effect)
  df2ggplot_all_genes_AUC <- data.frame()
  ## read AUC statistics files
  for (one_sample_size in sample_sizes) {
    for (one_replicate in 1:replicates) {
      dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, one_sample_size, one_sample_size, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
      ### Levene's test
      Levene_output_dir <- sprintf("%s/10-Levene/00-highly_DE/base_effect_%s/%s/repl_%d/00-all_genes", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      Levene_output_basename <- sprintf("%s_%s_Levene_test_results", dataset_name, normalization_method)
      auc_file <- sprintf("%s/%s_auc_stats.csv", Levene_output_dir, Levene_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("Levene", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
      ### MDSeq
      mdseq_output_dir <- sprintf("%s/20-MDSeq/00-highly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      mdseq_output_basename <- sprintf("%s_%s_MDSeq_results", dataset_name, normalization_method, DD_FC_threshold)
      auc_file <- sprintf("%s/%s_auc_stats.csv", mdseq_output_dir, mdseq_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("MDSeq", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
      ### DiPhiSeq
      diphiseq_output_dir <- sprintf("%s/30-DiPhiSeq/00-highly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      diphiseq_output_basename <- sprintf("%s_%s_DiPhiSeq_results", dataset_name, normalization_method)
      auc_file <- sprintf("%s/%s_auc_stats.csv", diphiseq_output_dir, diphiseq_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("DiPhiSeq", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
      ### GAMLSS
      gamlss_output_dir <- sprintf("%s/40-GAMLSS/00-highly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      gamlss_output_basename <- sprintf("%s_%s_GAMLSS_results", dataset_name, normalization_method)
      auc_file <- sprintf("%s/%s_auc_stats.csv", gamlss_output_dir, gamlss_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("GAMLSS", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
      ### DiffDist
      diffdist_output_dir <- sprintf("%s/50-DiffDist/00-highly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      diffdist_output_basename <- sprintf("%s_%s_DiffDist_results", dataset_name, normalization_method)
      auc_file <- sprintf("%s/%s_auc_stats.csv", diffdist_output_dir, diffdist_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("DiffDist", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
    }
  }
  if (nrow(df2ggplot_all_genes_AUC) > 0) {
    df2ggplot_all_genes_AUC$samp_cond_1 <- as.factor(df2ggplot_all_genes_AUC$samp_cond_1)
    ## order methods
    df2ggplot_all_genes_AUC$method <- factor(df2ggplot_all_genes_AUC$method, levels=method_order[method_order %in% df2ggplot_all_genes_AUC$method])
    ## change level names to be consistent with those of barplots
    levels(df2ggplot_all_genes_AUC$category) <- c("all genes", "highly DE", "lowly DE")
    
    ## plot
    plot_title <- sprintf("Differential dispersion performances - AUC - %d replicates \nDE: pct=%.1f, base effect=%.1f, effect param=%.2f, runif: %d;\nDV: pct=%.1f, base effect=%.1f, effect param=%.2f; runif: %d\nOutliers: single: high=%.2f, low=%.2f; random: high=%.2f, low=%.2f\nfilter: %d, outlier removal: %d\nMDSeq FC thresholds: mean=%.2f, dispersion=%.2f\nlabel FC threshold: mean=%.2f, dispersion=%.2f\np-val correction: DiPhiSeq:%s, MDSeq:%s", replicates, DE_fraction, DE_base_effect, DE_effect_param, non_DE_runif_FC, DD_fraction, DD_base_effect, DD_effect_param, non_DD_runif_FC, single_outlier_high_fraction, single_outlier_low_fraction, random_outlier_high_fraction, random_outlier_low_fraction, filter_threshold, outlier_removal_bool_param, DE_base_effect, DD_FC_threshold, DE_base_effect, DD_base_effect, diphiseq_pval_correction, mdseq_pval_correction)
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, plot_title, col = "black")
    p <- ggplot(data = df2ggplot_all_genes_AUC, aes(x=samp_cond_1, y=value, col=category)) + 
      geom_boxplot() +
      scale_color_brewer(palette="Dark2") +
      facet_wrap(~method, nrow=1) +
      theme_bw() +
      theme(panel.border = element_rect(color="grey50")) +
      theme(legend.position = "bottom") +
      theme(legend.title=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16)) +
      theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=14), strip.text = element_text(size=14)) +
      labs(x = "Sample size per condition", y = "AUC", fill="Differential mean\nexpression category")
    print(p)
  }
}
dev.off()
write.csv(df2ggplot_all_genes_AUC, file=sprintf("%s/%s.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)


# Figure 2: Ability to detect differential dispersion for lowly differentially expressed genes.
figure_basename <- "fig2"
## simulation parameters
### DE genes
DE_fraction <- 0
DE_effect_param <- 0

df2ggplot <- data.frame()
for (DE_base_effect in c(1.1, 1.2, 1.3, 1.4, 1.5)) {
  for (one_sample_size in sample_sizes) {
    for (one_replicate in 1:replicates) {
      dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, one_sample_size, one_sample_size, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
      ## read result performance files
      ### Levene's test
      Levene_output_dir <- sprintf("%s/10-Levene/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      Levene_output_basename <- sprintf("%s_%s_Levene_test_results", dataset_name, normalization_method)
      perf_file <- sprintf("%s/%s_performance_stats.csv", Levene_output_dir, Levene_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("Levene", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
      ### MDSeq
      mdseq_output_dir <- sprintf("%s/20-MDSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      mdseq_output_basename <- sprintf("%s_%s_MDSeq_results", dataset_name, normalization_method, DD_FC_threshold)
      perf_file <- sprintf("%s/%s_performance_stats.csv", mdseq_output_dir, mdseq_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("MDSeq", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
      ### DiPhiSeq
      diphiseq_output_dir <- sprintf("%s/30-DiPhiSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      diphiseq_output_basename <- sprintf("%s_%s_DiPhiSeq_results", dataset_name, normalization_method)
      perf_file <- sprintf("%s/%s_performance_stats.csv", diphiseq_output_dir, diphiseq_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("DiPhiSeq", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
      ### GAMLSS
      gamlss_output_dir <- sprintf("%s/40-GAMLSS/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      gamlss_output_basename <- sprintf("%s_%s_GAMLSS_results", dataset_name, normalization_method)
      perf_file <- sprintf("%s/%s_performance_stats.csv", gamlss_output_dir, gamlss_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("GAMLSS", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
      ### DiffDist
      diffdist_output_dir <- sprintf("%s/50-DiffDist/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      diffdist_output_basename <- sprintf("%s_%s_DiffDist_results", dataset_name, normalization_method)
      perf_file <- sprintf("%s/%s_performance_stats.csv", diffdist_output_dir, diffdist_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("DiffDist", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
    }
  }
}
pdf(sprintf("%s/%s.pdf", output_dir, figure_basename), width=9)
one_df2ggplot <- df2ggplot[which(df2ggplot$stat %in% c("FDR", "TPR")),]
one_df2ggplot$samp_cond_1 <- factor(one_df2ggplot$samp_cond_1)
one_df2ggplot$samp_cond_2 <- factor(one_df2ggplot$samp_cond_2)
one_df2ggplot$max_mean_FC <- factor(one_df2ggplot$max_mean_FC)
one_df2ggplot$method <- factor(one_df2ggplot$method, levels=method_order[method_order %in% one_df2ggplot$method])
p <- ggplot(data = one_df2ggplot, aes(x=samp_cond_1, y=value, col=max_mean_FC)) + 
  geom_boxplot() +
  scale_color_brewer(palette="Dark2") +
  facet_grid(stat ~ method, scales="free_y") +
  geom_hline(data=subset(one_df2ggplot, stat=="FDR"), aes(yintercept=0.05), linetype="dotted") +
  theme_bw() +
  theme(legend.position = "top") +
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16), legend.title=element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=12), strip.text = element_text(size=14)) +
  labs(x = "Sample size per condition", y="Value", col="Maximum mean fold-change value")
print(p)
dev.off()
write.csv(df2ggplot, file=sprintf("%s/%s.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)


# Figure 3: Differentially dispersed genes correctly identified by the evaluated methods among lowly differentially expressed genes.
figure_basename <- "fig3"
# path to output directory
output_dir <- "./output/simulations"
## simulation parameters
### DE genes
sample_sizes <- c(20, 30, 40, 50, 100)
DE_fraction <- 0
DE_effect_param <- 0
DE_base_effect <- 1.5
one_category <- "TP"
## plot parameters
correct_color <- brewer.pal(n = 8, name = "Dark2")[1]
opposite_color <- brewer.pal(n = 8, name = "Dark2")[2]
figure_barplot_theme <- theme(legend.direction="vertical") +
  theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=9)) +
  theme(axis.text.x=element_text(size=8, angle=20, hjust=1), axis.text.y=element_text(size=8), legend.text=element_text(size=7))
figure_scatterplot_theme <- theme_bw() +
  theme(panel.border = element_rect(color="grey50")) +
  theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=9)) +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=7))

logFC_consistency_same_signs <- function(results, results_col1, results_col2, annotations, annotations_col1, annotations_col2) {
  estimated_logFC <- results[[results_col1]] - results[[results_col2]]
  true_logFC <- log(annotations[rownames(results), annotations_col1] / annotations[rownames(results), annotations_col2], 2)
  logFC_consistency <- ifelse(true_logFC > 0 & estimated_logFC > 0 | true_logFC < 0 & estimated_logFC < 0, 1, 0)
  names(logFC_consistency) <- rownames(results)
  return(logFC_consistency)
}

logFC_consistency_opposite_signs <- function(results, results_col, annotations, annotations_col1, annotations_col2) {
  estimated_logFC <- results[[results_col]]
  true_logFC <- log(annotations[rownames(results), annotations_col1] / annotations[rownames(results), annotations_col2], 2)
  logFC_consistency <- ifelse(true_logFC > 0 & estimated_logFC < 0 | true_logFC < 0 & estimated_logFC > 0, 1, 0)
  names(logFC_consistency) <- rownames(results)
  return(logFC_consistency)
}

mean_counts_over_methods_categories <- function(counts) {
  mean_counts <- counts %>%
    group_by(across(all_of(c("method", "category")))) %>%
    summarise(mean_count = mean(count, na.rm=TRUE))
  # add mean percentages
  pct_vector <- c()
  mean_counts <- as.data.frame(mean_counts)
  for (one_method in levels(mean_counts$method)) {
    pct_vector <- c(pct_vector, round(mean_counts[which(mean_counts$method==one_method), "mean_count"] / sum(mean_counts[which(mean_counts$method==one_method), "mean_count"], na.rm=TRUE)*100, 2))
  }
  mean_counts <- cbind(mean_counts, pct=pct_vector)
  return(mean_counts)
}


all_sample_sizes_mean_counts_per_method_df <- data.frame()
all_sample_sizes_all_replicates_logFC_mean_counts_df <- data.frame()
all_sample_sizes_all_replicates_wrong_diffdist_logFC_mean_counts_df <- data.frame()
all_sample_sizes_all_replicates_real_DiffDist_GAMLSS_logFC_df <- data.frame()
pdf(sprintf("%s/fig3.pdf", output_dir), width=9)
for (nb_samples in sample_sizes) {
  print(sprintf("nb samples: %d", nb_samples))
  nb_samples_condition_1 <- nb_samples
  nb_samples_condition_2 <- nb_samples
  
  # data frames
  ## upset
  all_replicates_genes_methods_df <- data.frame()
  all_replicate_genes_method_list_df <- data.frame()
  all_replicates_counts_per_method_df <- data.frame()
  ## logFC consistency
  all_replicates_logFC_counts_df <- data.frame()
  all_replicates_wrong_diffdist_logFC_counts_df <- data.frame()
  ## real and DiffDist dispersion log FC
  all_replicates_real_DiffDist_logFC_df <- data.frame()
  ## real and GAMLSS dispersion log FC
  all_replicates_real_GAMLSS_logFC_df  <- data.frame()
  
  for (one_replicate in seq(1:replicates)) {
    ## read result performance files
    dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, nb_samples, nb_samples, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
    ### Levene's test
    Levene_output_dir <- sprintf("%s/10-Levene/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), nb_samples, one_replicate)
    Levene_results_file <- sprintf("%s/%s_%s_Levene_test_results_performance.csv", Levene_output_dir, dataset_name, normalization_method)
    if (file.exists(Levene_results_file)) {
      Levene_results_df <- read.csv(Levene_results_file, header=TRUE, row.names=1)
    }
    ### MDSeq
    mdseq_output_dir <- sprintf("%s/20-MDSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), nb_samples, one_replicate)
    mdseq_results_file <- sprintf("%s/%s_%s_MDSeq_results_performance.csv", mdseq_output_dir, dataset_name, normalization_method, DD_FC_threshold)
    if (file.exists(mdseq_results_file)) {
      mdseq_results_df <- read.csv(mdseq_results_file, header=TRUE, row.names=1)
    }
    ### DiPhiSeq
    diphiseq_output_dir <- sprintf("%s/30-DiPhiSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), nb_samples, one_replicate)
    diphiseq_results_file <- sprintf("%s/%s_%s_DiPhiSeq_results_performance.csv", diphiseq_output_dir, dataset_name, normalization_method)
    if (file.exists(diphiseq_results_file)) {
      diphiseq_results_df <- read.csv(diphiseq_results_file, header=TRUE, row.names=1)
    }
    ### GAMLSS
    gamlss_output_dir <- sprintf("%s/40-GAMLSS/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), nb_samples, one_replicate)
    gamlss_results_file <- sprintf("%s/%s_%s_GAMLSS_results_performance.csv", gamlss_output_dir, dataset_name, normalization_method)
    if (file.exists(gamlss_results_file)) {
      gamlss_results_df <- read.csv(gamlss_results_file, header=TRUE, row.names=1)
    }
    ### DiffDist
    diffdist_output_dir <- sprintf("%s/50-DiffDist/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), nb_samples, one_replicate)
    diffdist_results_file <- sprintf("%s/%s_%s_DiffDist_results_performance.csv", diffdist_output_dir, dataset_name, normalization_method)
    if (file.exists(diffdist_results_file)) {
      diffdist_results_df <- read.csv(diffdist_results_file, header=TRUE, row.names=1)
    }
    
    # get TP genes
    one_replicate_genes_methods_df <- data.frame()
    if (file.exists(Levene_results_file)) {
      one_replicate_genes_methods_df <- rbind(one_replicate_genes_methods_df, data.frame(gene=rownames(Levene_results_df[which(Levene_results_df$category.p_value.BH==one_category),]), method=rep("Levene", length(rownames(Levene_results_df[which(Levene_results_df$category.p_value.BH==one_category),])))))
    }
    if (file.exists(mdseq_results_file)) {
      one_replicate_genes_methods_df <- rbind(one_replicate_genes_methods_df, data.frame(gene=rownames(mdseq_results_df[which(mdseq_results_df$category.FDR.dispersion==one_category),]), method=rep("MDSeq", length(rownames(mdseq_results_df[which(mdseq_results_df$category.FDR.dispersion==one_category),])))))
    }
    if (file.exists(diphiseq_results_file)) {
      one_replicate_genes_methods_df <- rbind(one_replicate_genes_methods_df, data.frame(gene=rownames(diphiseq_results_df[which(diphiseq_results_df$category.fdr.phi==one_category),]), method=rep("DiPhiSeq", length(rownames(diphiseq_results_df[which(diphiseq_results_df$category.fdr.phi==one_category),])))))
    }
    if (file.exists(gamlss_results_file)) {
      one_replicate_genes_methods_df <- rbind(one_replicate_genes_methods_df, data.frame(gene=rownames(gamlss_results_df[which(gamlss_results_df$category.padj.cv==one_category),]), method=rep("GAMLSS", length(rownames(gamlss_results_df[which(gamlss_results_df$category.padj.cv==one_category),])))))
    }
    
    if (file.exists(diffdist_results_file)) {
      one_replicate_genes_methods_df <- rbind(one_replicate_genes_methods_df, data.frame(gene=rownames(diffdist_results_df[which(diffdist_results_df$category.disp.1vs2.pval.BH==one_category),]), method=rep("DiffDist", length(rownames(diffdist_results_df[which(diffdist_results_df$category.disp.1vs2.pval.BH==one_category),])))))
    }
    
    if (dim(one_replicate_genes_methods_df)[1] > 0) {
      one_replicate_genes_methods_df$gene <- as.character(one_replicate_genes_methods_df$gene)
      one_replicate_genes_methods_df$method <- as.character(one_replicate_genes_methods_df$method)
      one_replicate_genes_methods_df <- cbind(replicate=rep(sprintf("repl_%d", one_replicate), dim(one_replicate_genes_methods_df)[1]), one_replicate_genes_methods_df)
      all_replicates_genes_methods_df <- rbind(all_replicates_genes_methods_df, one_replicate_genes_methods_df)
      # list of methods per gene
      one_replicate_genes_method_list_df <- one_replicate_genes_methods_df %>%
        group_by(gene) %>%
        summarise(method = list(method))
      one_replicate_genes_method_list_df <- one_replicate_genes_method_list_df %>%
        group_by(gene) %>%
        summarise(method = unlist(lapply(method, function(x) {paste(x, collapse="&")})))
      all_replicate_genes_method_list_df <- rbind(all_replicate_genes_method_list_df, one_replicate_genes_method_list_df)
      # counts per method
      counts_per_method <- table(one_replicate_genes_method_list_df$method)
      one_replicate_counts_per_method_df <- data.frame(replicate=rep(sprintf("repl_%d", one_replicate), length(counts_per_method)), method=names(counts_per_method), count=as.numeric(counts_per_method))
      all_replicates_counts_per_method_df <- rbind(all_replicates_counts_per_method_df, one_replicate_counts_per_method_df)
      
      # load simulated data to get true mean and dispersion fold-changes
      data_dir <- sprintf("%s/00-Data/10-lowly_DE/base_effect_%s/%d/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), nb_samples, one_replicate)
      # dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, nb_samples, nb_samples, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
      rds_file <- sprintf("%s/%s.rds", data_dir, dataset_name)
      if (file.exists(rds_file)) {
        sim_data <- readRDS(rds_file) # contains compData object (compcodeR)
        ## simulation parameters
        conditions <- as.factor(sim_data@sample.annotations$condition)
        annotations_df <- sim_data@variable.annotations
        
        # logFC consistency
        ## GAMLSS
        if (file.exists(gamlss_results_file)) {
          gamlss_TP_results_df <- gamlss_results_df[which(gamlss_results_df$category.padj.cv==sprintf("%s", one_category)),]
          TP_true_mean_FC <- annotations_df[rownames(gamlss_TP_results_df), "truelog2foldchanges"]
          TP_true_disp_FC <- log(annotations_df[rownames(gamlss_TP_results_df), "truedispersions.S2"] / annotations_df[rownames(gamlss_TP_results_df), "truedispersions.S1"], 2)
          gamlss_TP_results_df <- cbind(gamlss_TP_results_df, trueFC.mean=TP_true_mean_FC, trueFC.disp=TP_true_disp_FC)
          gamlss_TP_logFC_consistency <- logFC_consistency_same_signs(gamlss_TP_results_df, "CV.2", "CV.1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          if ("0" %in% names(table(gamlss_TP_logFC_consistency))) {
            gamlss_TP_logFC_consistency_table_counts <- as.numeric(unname(table(gamlss_TP_logFC_consistency)))
          } else {
            gamlss_TP_logFC_consistency_table_counts <- c(0, as.numeric(unname(table(gamlss_TP_logFC_consistency))))
          }
          ### get genes with wrong logFCs
          gamlss_TP_results_df <- cbind(gamlss_TP_results_df, dispersion_logFC_consistency=gamlss_TP_logFC_consistency)
          wrong_gamlss_TP_logFC_genes <- rownames(gamlss_TP_results_df[which(gamlss_TP_results_df$dispersion_logFC_consistency==0),])
        } else {
          gamlss_TP_logFC_consistency_table_counts <- c(NA, NA)
          wrong_gamlss_TP_logFC_genes <- NA
        }
        
        ## DiffDist
        if (file.exists(diffdist_results_file)) {
          diffdist_TP_results_df <- diffdist_results_df[which(diffdist_results_df$category.disp.1vs2.pval.BH==sprintf("%s", one_category)),]
          diffdist_TP_logFC_consistency <- logFC_consistency_opposite_signs(diffdist_TP_results_df, sprintf("disp.%svs%s.logFC", levels(conditions)[1], levels(conditions)[2]), annotations_df, "truedispersions.S2", "truedispersions.S1")
          if ("0" %in% names(table(diffdist_TP_logFC_consistency))) {
            diffdist_TP_logFC_consistency_table_counts <- as.numeric(unname(table(diffdist_TP_logFC_consistency)))
          } else {
            diffdist_TP_logFC_consistency_table_counts <- c(0, as.numeric(unname(table(diffdist_TP_logFC_consistency))))
          }
          ### get genes with wrong logFCs
          diffdist_TP_results_df <- cbind(diffdist_TP_results_df, dispersion_logFC_consistency=diffdist_TP_logFC_consistency)
          wrong_diffdist_TP_logFC_genes <- rownames(diffdist_TP_results_df[which(diffdist_TP_results_df$dispersion_logFC_consistency==0),])
        } else {
          diffdist_TP_logFC_consistency_table_counts <- c(NA, NA)
          wrong_diffdist_TP_logFC_genes <- NA
        }
        
        ## Levene test
        if (file.exists(Levene_results_file)) {
          Levene_test_TP_results_df <- Levene_results_df[which(Levene_results_df$category.p_value.BH==sprintf("%s", one_category)),]
          Levene_test_TP_logFC_consistency <- logFC_consistency_same_signs(Levene_test_TP_results_df, "var_2", "var_1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          if ("0" %in% names(table(Levene_test_TP_logFC_consistency))) {
            Levene_test_TP_logFC_consistency_table_counts <- as.numeric(unname(table(Levene_test_TP_logFC_consistency)))
          } else {
            Levene_test_TP_logFC_consistency_table_counts <- c(0, as.numeric(unname(table(Levene_test_TP_logFC_consistency))))
          }
          ### get log FC consistency for DiffDist and GAMLSS genes with wrong logFCs
          if (! is.na(wrong_gamlss_TP_logFC_genes)) {
            Levene_test_wrong_gamlss_logFC_consistency <- logFC_consistency_same_signs(Levene_results_df[wrong_gamlss_TP_logFC_genes,], "var_2", "var_1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          } else {
            Levene_test_wrong_gamlss_logFC_consistency <- c(NA, NA)
          }
          if (! is.na(wrong_diffdist_TP_logFC_genes)) {
            Levene_test_wrong_diffdist_logFC_consistency <- logFC_consistency_same_signs(Levene_results_df[wrong_diffdist_TP_logFC_genes,], "var_2", "var_1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          } else {
            Levene_test_wrong_diffdist_logFC_consistency <- c(NA, NA)
          }
        } else {
          Levene_test_TP_logFC_consistency_table_counts <- c(NA, NA)
          Levene_test_wrong_gamlss_logFC_consistency <- c(NA, NA)
          Levene_test_wrong_diffdist_logFC_consistency <- c(NA, NA)
        }
        
        ## MDSeq
        if (file.exists(mdseq_results_file)) {
          mdseq_TP_results_df <- mdseq_results_df[which(mdseq_results_df$category.FDR.dispersion==sprintf("%s", one_category)),]
          mdseq_TP_logFC_consistency <- logFC_consistency_opposite_signs(mdseq_TP_results_df, sprintf("X%svs%s.dispersion.log2FC.0", levels(conditions)[1], levels(conditions)[2]), annotations_df, "truedispersions.S2", "truedispersions.S1")
          if ("0" %in% names(table(mdseq_TP_logFC_consistency))) {
            mdseq_TP_logFC_consistency_table_counts <- as.numeric(unname(table(mdseq_TP_logFC_consistency)))
          } else {
            mdseq_TP_logFC_consistency_table_counts <- c(0, as.numeric(unname(table(mdseq_TP_logFC_consistency))))
          }
          ### get log FC consistency for DiffDist and GAMLSS genes with wrong logFCs
          if (! is.na(wrong_gamlss_TP_logFC_genes)) {
            mdseq_wrong_gamlss_logFC_consistency <- logFC_consistency_opposite_signs(mdseq_results_df[wrong_gamlss_TP_logFC_genes,], sprintf("X%svs%s.dispersion.log2FC.0", levels(conditions)[1], levels(conditions)[2]), annotations_df, "truedispersions.S2", "truedispersions.S1")
          } else {
            mdseq_wrong_gamlss_logFC_consistency <- c(NA, NA)
          }
          if (! is.na(wrong_diffdist_TP_logFC_genes)) {
            mdseq_wrong_diffdist_logFC_consistency <- logFC_consistency_opposite_signs(mdseq_results_df[wrong_diffdist_TP_logFC_genes,], sprintf("X%svs%s.dispersion.log2FC.0", levels(conditions)[1], levels(conditions)[2]), annotations_df, "truedispersions.S2", "truedispersions.S1")
          } else {
            mdseq_wrong_diffdist_logFC_consistency <- c(NA, NA)
          }
        } else {
          mdseq_TP_logFC_consistency_table_counts <- c(NA, NA)
          mdseq_wrong_gamlss_logFC_consistency <- c(NA, NA)
          mdseq_wrong_diffdist_logFC_consistency <- c(NA, NA)
        }
        
        ## DiPhiSeq
        if (file.exists(diphiseq_results_file)) {
          diphiseq_TP_results_df <- diphiseq_results_df[which(diphiseq_results_df$category.fdr.phi==sprintf("%s", one_category)),]
          diphiseq_TP_logFC_consistency <- logFC_consistency_same_signs(diphiseq_TP_results_df, "phi2", "phi1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          if ("0" %in% names(table(diphiseq_TP_logFC_consistency))) {
            diphiseq_TP_logFC_consistency_table_counts <- as.numeric(unname(table(diphiseq_TP_logFC_consistency)))
          } else {
            if ("1" %in% names(table(diphiseq_TP_logFC_consistency))) {
              diphiseq_TP_logFC_consistency_table_counts <- c(0, as.numeric(unname(table(diphiseq_TP_logFC_consistency))))
            } else {
              diphiseq_TP_logFC_consistency_table_counts <- c(0, 0)
            }
          }
          ### get log FC consistency for DiffDist and GAMLSS genes with wrong logFCs
          if (! is.na(wrong_gamlss_TP_logFC_genes)) {
            diphiseq_wrong_gamlss_logFC_consistency <- logFC_consistency_same_signs(diphiseq_results_df[wrong_gamlss_TP_logFC_genes,], "phi2", "phi1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          } else {
            diphiseq_wrong_gamlss_logFC_consistency <- c(NA, NA)
          }
          if (! is.na(wrong_diffdist_TP_logFC_genes)) {
            diphiseq_wrong_diffdist_logFC_consistency <- logFC_consistency_same_signs(diphiseq_results_df[wrong_diffdist_TP_logFC_genes,], "phi2", "phi1", annotations_df, "truedispersions.S2", "truedispersions.S1")
          } else {
            diphiseq_wrong_diffdist_logFC_consistency <- c(NA, NA)
          }
        } else {
          diphiseq_TP_logFC_consistency_table_counts <- c(NA, NA)
          diphiseq_wrong_gamlss_logFC_consistency <- c(NA, NA)
          diphiseq_wrong_diffdist_logFC_consistency <- c(NA, NA)
        }
        
        one_replicate_counts_df <- data.frame(replicate=rep(one_replicate, 10), 
                                              method=c(rep("DiffDist", 2), rep("DiPhiSeq", 2), rep("MDSeq", 2), rep("Levene", 2), rep("GAMLSS", 2)), 
                                              category=rep(c("opposite", "correct"), 5), 
                                              count=c(diffdist_TP_logFC_consistency_table_counts, diphiseq_TP_logFC_consistency_table_counts, mdseq_TP_logFC_consistency_table_counts, Levene_test_TP_logFC_consistency_table_counts, gamlss_TP_logFC_consistency_table_counts), 
                                              pct=c(round(diffdist_TP_logFC_consistency_table_counts/sum(diffdist_TP_logFC_consistency_table_counts)*100, 2), round(diphiseq_TP_logFC_consistency_table_counts/sum(diphiseq_TP_logFC_consistency_table_counts)*100, 2), round(mdseq_TP_logFC_consistency_table_counts/sum(mdseq_TP_logFC_consistency_table_counts)*100, 2), round(Levene_test_TP_logFC_consistency_table_counts/sum(Levene_test_TP_logFC_consistency_table_counts)*100, 2), round(gamlss_TP_logFC_consistency_table_counts/sum(gamlss_TP_logFC_consistency_table_counts)*100, 2)))
        all_replicates_logFC_counts_df <- rbind(all_replicates_logFC_counts_df, one_replicate_counts_df)
        
        one_replicate_wrong_diffdist_logFC_counts_df <- data.frame(replicate=rep(one_replicate, 12), 
                                                                   wrong=c(rep("DiffDist", 6), rep("GAMLSS", 6)),
                                                                   method=rep(c(rep("DiPhiSeq", 2), rep("MDSeq", 2), rep("Levene", 2)), 2),
                                                                   category=rep(c("opposite", "correct"), 6), 
                                                                   count=c(table(diphiseq_wrong_diffdist_logFC_consistency)["0"], table(diphiseq_wrong_diffdist_logFC_consistency)["1"], table(mdseq_wrong_diffdist_logFC_consistency)["0"], table(mdseq_wrong_diffdist_logFC_consistency)["1"], table(Levene_test_wrong_diffdist_logFC_consistency)["0"], table(Levene_test_wrong_diffdist_logFC_consistency)["1"],
                                                                           table(diphiseq_wrong_gamlss_logFC_consistency)["0"], table(diphiseq_wrong_gamlss_logFC_consistency)["1"], table(mdseq_wrong_gamlss_logFC_consistency)["0"], table(mdseq_wrong_gamlss_logFC_consistency)["1"], table(Levene_test_wrong_gamlss_logFC_consistency)["0"], table(Levene_test_wrong_gamlss_logFC_consistency)["1"]), 
                                                                   pct=c(round(c(table(diphiseq_wrong_diffdist_logFC_consistency)["0"], table(diphiseq_wrong_diffdist_logFC_consistency)["1"])/sum(table(diphiseq_wrong_diffdist_logFC_consistency))*100, 2), round(c(table(mdseq_wrong_diffdist_logFC_consistency)["0"], table(mdseq_wrong_diffdist_logFC_consistency)["1"])/sum(table(mdseq_wrong_diffdist_logFC_consistency))*100, 2), round(c(table(Levene_test_wrong_diffdist_logFC_consistency)["0"], table(Levene_test_wrong_diffdist_logFC_consistency)["1"])/sum(table(Levene_test_wrong_diffdist_logFC_consistency))*100, 2),
                                                                         round(c(table(diphiseq_wrong_gamlss_logFC_consistency)["0"], table(diphiseq_wrong_gamlss_logFC_consistency)["1"])/sum(table(diphiseq_wrong_gamlss_logFC_consistency))*100, 2), round(c(table(mdseq_wrong_gamlss_logFC_consistency)["0"], table(mdseq_wrong_gamlss_logFC_consistency)["1"])/sum(table(mdseq_wrong_gamlss_logFC_consistency))*100, 2), round(c(table(Levene_test_wrong_gamlss_logFC_consistency)["0"], table(Levene_test_wrong_gamlss_logFC_consistency)["1"])/sum(table(Levene_test_wrong_gamlss_logFC_consistency))*100, 2)))
        all_replicates_wrong_diffdist_logFC_counts_df <- rbind(all_replicates_wrong_diffdist_logFC_counts_df, one_replicate_wrong_diffdist_logFC_counts_df)
        
        # real and GAMLSS dispersion logFC scatter plots
        if (file.exists(gamlss_results_file)) {
          true_mean_fcs <- log(annotations_df[rownames(gamlss_TP_results_df), "truemeans.S2"] / annotations_df[rownames(gamlss_TP_results_df), "truemeans.S1"], 2)
          true_dispersion_fcs <- log(annotations_df[rownames(gamlss_TP_results_df), "truedispersions.S2"] / annotations_df[rownames(gamlss_TP_results_df), "truedispersions.S1"], 2)
          
          df2ggplot <- data.frame(replicate=rep(one_replicate, dim(gamlss_TP_results_df)[1]), est_FC_mean=log(gamlss_TP_results_df[["cpm.2"]]/gamlss_TP_results_df[["cpm.1"]], 2), est_FC_dispersion=log(gamlss_TP_results_df[["CV.2"]]/gamlss_TP_results_df[["CV.1"]], 2), true_FC_mean=true_mean_fcs, true_FC_dispersion=true_dispersion_fcs, category.logFC.dispersion_consistency=gamlss_TP_results_df[["dispersion_logFC_consistency"]])
          df2ggplot$category.logFC.dispersion_consistency <- as.factor(df2ggplot$category.logFC.dispersion_consistency)
          all_replicates_real_GAMLSS_logFC_df <- rbind(all_replicates_real_GAMLSS_logFC_df, df2ggplot)
        }
        
        # real and DiffDist dispersion logFC scatter plots
        if (file.exists(diffdist_results_file)) {
          true_mean_fcs <- log(annotations_df[rownames(diffdist_TP_results_df), "truemeans.S2"] / annotations_df[rownames(diffdist_TP_results_df), "truemeans.S1"], 2)
          true_dispersion_fcs <- log(annotations_df[rownames(diffdist_TP_results_df), "truedispersions.S2"] / annotations_df[rownames(diffdist_TP_results_df), "truedispersions.S1"], 2)
          df2ggplot <- data.frame(replicate=rep(one_replicate, dim(diffdist_TP_results_df)[1]), est_FC_mean=-diffdist_TP_results_df[["mean.1vs2.logFC"]], est_FC_dispersion=-diffdist_TP_results_df[["disp.1vs2.logFC"]], true_FC_mean=true_mean_fcs, true_FC_dispersion=true_dispersion_fcs, category.logFC.dispersion_consistency=diffdist_TP_results_df[["dispersion_logFC_consistency"]])
          df2ggplot$category.logFC.dispersion_consistency <- as.factor(df2ggplot$category.logFC.dispersion_consistency)
          all_replicates_real_DiffDist_logFC_df <- rbind(all_replicates_real_DiffDist_logFC_df, df2ggplot)
        }
      }
    }
  }
  
  if (dim(all_replicates_counts_per_method_df)[1] > 0) {
    # write outputs
    output_basename <- sprintf("DE_%s_0_runif_%d_%d_%s", gsub("[.]", "_", DE_base_effect), nb_samples_condition_1, nb_samples_condition_2, one_category)
    
    # mean counts plots
    ## upset
    mean_counts_per_method_df <- all_replicates_counts_per_method_df %>%
      group_by(method) %>%
      summarise(mean_count = mean(count))
    all_sample_sizes_mean_counts_per_method_df <- rbind(all_sample_sizes_mean_counts_per_method_df, cbind(sample_size=rep(nb_samples, dim(mean_counts_per_method_df)[1]), mean_counts_per_method_df))
    expression_upset_input <- mean_counts_per_method_df$mean_count
    names(expression_upset_input) <- mean_counts_per_method_df$method
    upset_plot_figure <- upset(fromExpression(expression_upset_input), order.by = "freq", point.size=2, line.size=1, mb.ratio=c(0.6, 0.4))
    
    ## logFC consistency
    ### all TP results
    all_replicates_logFC_mean_counts_df <- mean_counts_over_methods_categories(all_replicates_logFC_counts_df)
    all_sample_sizes_all_replicates_logFC_mean_counts_df <- rbind(all_sample_sizes_all_replicates_logFC_mean_counts_df, cbind(sample_size=rep(nb_samples, dim(all_replicates_logFC_mean_counts_df)[1]), all_replicates_logFC_mean_counts_df))
    ### DiffDist and GAMLSS TP genes with wrong logFCs
    all_replicates_wrong_diffdist_logFC_mean_counts_df <- data.frame()
    for (one_method in levels(all_replicates_wrong_diffdist_logFC_counts_df$wrong)) {
      one_method_mean_counts_df <- mean_counts_over_methods_categories(all_replicates_wrong_diffdist_logFC_counts_df[which(all_replicates_wrong_diffdist_logFC_counts_df$wrong==one_method),])
      all_replicates_wrong_diffdist_logFC_mean_counts_df <- rbind(all_replicates_wrong_diffdist_logFC_mean_counts_df, cbind(wrong=rep(one_method, dim(one_method_mean_counts_df)[1]), one_method_mean_counts_df))
    }
    all_sample_sizes_all_replicates_wrong_diffdist_logFC_mean_counts_df <- rbind(all_sample_sizes_all_replicates_wrong_diffdist_logFC_mean_counts_df, cbind(sample_size=rep(nb_samples, dim(all_replicates_wrong_diffdist_logFC_mean_counts_df)[1]), all_replicates_wrong_diffdist_logFC_mean_counts_df))
    all_replicates_logFC_mean_counts_df$method <- factor(all_replicates_logFC_mean_counts_df$method, levels=method_order[method_order %in% all_replicates_logFC_mean_counts_df$method])
    logFC_consistency_plot_figure <- ggplot(all_replicates_logFC_mean_counts_df, aes(x=method, y=mean_count, fill=category)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=pct), size=3, position=position_stack(vjust=0.5)) +
      labs(x="Method", y="#genes", fill="Dispersion\nlog fold-change sign") +
      scale_fill_manual(values=c(correct_color, opposite_color)) +
      theme_bw() +
      theme(legend.position = "bottom") +
      theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10)) 
    all_replicates_wrong_diffdist_logFC_mean_counts_df$method <- factor(all_replicates_wrong_diffdist_logFC_mean_counts_df$method, levels=method_order[method_order %in% all_replicates_wrong_diffdist_logFC_mean_counts_df$method])
    all_replicates_wrong_diffdist_logFC_mean_counts_df$wrong <- factor(all_replicates_wrong_diffdist_logFC_mean_counts_df$wrong, levels=method_order[method_order %in% all_replicates_wrong_diffdist_logFC_mean_counts_df$wrong])
    
    if (sum(is.na(all_replicates_wrong_diffdist_logFC_mean_counts_df$mean_count)) == length(all_replicates_wrong_diffdist_logFC_mean_counts_df$mean_count)) {
      # empty plot
      wrong_diffdist_logFC_consistency_plot_figure <- ggplot() + theme_void()
    } else {
      wrong_diffdist_logFC_consistency_plot_figure <- ggplot(all_replicates_wrong_diffdist_logFC_mean_counts_df, aes(x=method, y=mean_count, fill=category)) +
        geom_bar(stat="identity") +
        facet_wrap(~wrong, nrow=2) +
        geom_text(aes(label=pct), size=3, position=position_stack(vjust=0.5)) +
        labs(x="Method", y="#genes", fill="Dispersion\nlog fold-change sign") +
        scale_fill_manual(values=c(correct_color, opposite_color)) +
        theme_bw() +
        theme(legend.position = "bottom") +
        theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12)) +
        theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10))
    }
    
    if (dim(all_replicates_real_DiffDist_logFC_df)[1] > 0) {
      ## change level names to be consistent with those of barplots
      levels(all_replicates_real_DiffDist_logFC_df$category.logFC.dispersion_consistency) <- c("opposite", "correct")
      ## change level order
      all_replicates_real_DiffDist_logFC_df$category.logFC.dispersion_consistency <- factor(all_replicates_real_DiffDist_logFC_df$category.logFC.dispersion_consistency, levels=c("correct", "opposite"))
    }
    if (dim(all_replicates_real_GAMLSS_logFC_df)[1] > 0) {
      ## change level names to be consistent with those of barplots
      levels(all_replicates_real_GAMLSS_logFC_df$category.logFC.dispersion_consistency) <- c("opposite", "correct")
      ## change level order
      all_replicates_real_GAMLSS_logFC_df$category.logFC.dispersion_consistency <- factor(all_replicates_real_GAMLSS_logFC_df$category.logFC.dispersion_consistency, levels=c("correct", "opposite"))
      ### remove extreme values for plot
      all_replicates_real_GAMLSS_logFC_df2plot <- all_replicates_real_GAMLSS_logFC_df[which(all_replicates_real_GAMLSS_logFC_df$est_FC_dispersion > -5 & all_replicates_real_GAMLSS_logFC_df$est_FC_dispersion < 5),]
    } else {
      # empty data frame
      all_replicates_real_GAMLSS_logFC_df2plot <- data.frame()
    }
    all_replicates_real_DiffDist_GAMLSS_logFC_df <- rbind(cbind(method=rep("DiffDist", dim(all_replicates_real_DiffDist_logFC_df)[1]), all_replicates_real_DiffDist_logFC_df), cbind(method=rep("GAMLSS", dim(all_replicates_real_GAMLSS_logFC_df2plot)[1]), all_replicates_real_GAMLSS_logFC_df2plot))
    if (dim(all_replicates_real_DiffDist_GAMLSS_logFC_df)[1] > 0) {
      all_sample_sizes_all_replicates_real_DiffDist_GAMLSS_logFC_df <- rbind(all_sample_sizes_all_replicates_real_DiffDist_GAMLSS_logFC_df, cbind(sample_size=rep(nb_samples, dim(all_replicates_real_DiffDist_GAMLSS_logFC_df)[1]), all_replicates_real_DiffDist_GAMLSS_logFC_df))
      
      all_replicates_real_DiffDist_GAMLSS_logFC_df$method <- factor(all_replicates_real_DiffDist_GAMLSS_logFC_df$method, levels=method_order[method_order %in% all_replicates_real_DiffDist_GAMLSS_logFC_df$method])
      
      true_fcs_figure <- ggplot(all_replicates_real_DiffDist_GAMLSS_logFC_df, aes(x=true_FC_mean, y=true_FC_dispersion, col=category.logFC.dispersion_consistency)) +
        geom_point(size=0.5) +
        geom_point(data=subset(all_replicates_real_DiffDist_logFC_df, category.logFC.dispersion_consistency=="opposite"), aes(x=true_FC_mean, y=true_FC_dispersion, col=category.logFC.dispersion_consistency), size=0.5) + # faire resortir les "0"
        facet_wrap(~method, nrow=2, scales="free_y") +
        scale_color_manual(values=c(correct_color, opposite_color)) +
        geom_hline(yintercept=0, linetype="solid", color="black", alpha=0.3) +
        geom_vline(xintercept=0, linetype="solid", color="black", alpha=0.3) +
        labs(x="Real mean fold-change", y="Real dispersion fold-change", col="Dispersion log fold-change sign") +
        figure_scatterplot_theme +
        theme(legend.position = "bottom")
      
      disp_fcs_figure <- ggplot(all_replicates_real_DiffDist_GAMLSS_logFC_df, aes(x=true_FC_dispersion, y=est_FC_dispersion, col=category.logFC.dispersion_consistency)) +
        geom_point(size=0.5) +
        facet_wrap(~method, nrow=2, scales="free_y") +
        scale_color_manual(values=c(correct_color, opposite_color)) +
        geom_hline(yintercept=0, linetype="solid", color="black", alpha=0.3) +
        geom_vline(xintercept=0, linetype="solid", color="black", alpha=0.3) +
        labs(x="Real dispersion fold-change", y="Estimated dispersion fold-change", col="Dispersion log fold-change sign") +
        figure_scatterplot_theme + 
        theme(legend.position = "bottom")
    } else {
      # empty plots
      true_fcs_figure <- disp_fcs_figure <- ggplot() + theme_void()
    }
    
    # convert UpSetR plot to ggplot plot: https://github.com/hms-dbmi/UpSetR/issues/105
    upset_plot_figure_ggplot <- cowplot::plot_grid(NULL, upset_plot_figure$Main_bar, upset_plot_figure$Sizes, upset_plot_figure$Matrix,
                                                   nrow=2, align='hv', rel_heights = c(2.5,1),
                                                   rel_widths = c(1,3))
    logFC_consistency_plot_figure <- logFC_consistency_plot_figure + figure_barplot_theme
    wrong_diffdist_logFC_consistency_plot_figure <- wrong_diffdist_logFC_consistency_plot_figure + figure_barplot_theme
    
    plot_legend_var <- get_legend(true_fcs_figure)
    true_fcs_figure_no_legend <- true_fcs_figure + theme(legend.position = "none")
    disp_fcs_figure_no_legend <- disp_fcs_figure + theme(legend.position = "none")
    real_DiffDist_GAMLSS_logFC_plot_figure <- ggdraw() +
      draw_plot(true_fcs_figure_no_legend, 0, 0.1, 0.5, 0.9) +
      draw_plot(disp_fcs_figure_no_legend, 0.5, 0.1, 0.5, 0.9) +
      draw_plot(plot_legend_var, 0, 0, 1, 0.1)
    
    figure_plot <- ggdraw() +
      draw_plot(upset_plot_figure_ggplot, 0, 0.5, 1, 0.5) +
      draw_plot(logFC_consistency_plot_figure, 0, 0, 0.25, 0.5) +
      draw_plot(real_DiffDist_GAMLSS_logFC_plot_figure, 0.25, 0, 0.5, 0.5) +
      # draw_plot(disp_fcs_figure, 0.495, 0, 0.33, 0.5) +
      draw_plot(wrong_diffdist_logFC_consistency_plot_figure, 0.75, 0, 0.25, 0.5) +
      draw_plot_label(c("A", "B", "C", "D"), c(0, 0, 0.25, 0.75), c(1, 0.5, 0.5, 0.5), size=15)
    
    plot_title <- ggdraw() + draw_label(sprintf("Sample size: %d", nb_samples), fontface="bold", size=16)
    figure_plot_title <- plot_grid(plot_title, figure_plot, ncol=1, rel_heights=c(0.03, 0.97))
    print(figure_plot_title)
  }
}
dev.off()
write.csv(all_sample_sizes_mean_counts_per_method_df, file=sprintf("%s/fig3_upset_methods_mean_counts.csv", output_dir), quote=FALSE, row.names=FALSE)
write.csv(all_sample_sizes_all_replicates_logFC_mean_counts_df, file=sprintf("%s/fig3_logFC_signs_mean_counts.csv", output_dir), quote=FALSE, row.names=FALSE)
write.csv(all_sample_sizes_all_replicates_wrong_diffdist_logFC_mean_counts_df, file=sprintf("%s/fig3_wrong_diffdist_logFC_mean_counts.csv", output_dir), quote=FALSE, row.names=FALSE)
write.csv(all_sample_sizes_all_replicates_real_DiffDist_GAMLSS_logFC_df, file=sprintf("%s/fig3_real_DiffDist_GAMLSS_dispersion_logFC.csv", output_dir), quote=FALSE, row.names=FALSE)


