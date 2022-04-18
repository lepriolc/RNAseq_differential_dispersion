library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(VennDiagram)


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

# Figure 1: DiPhiSeq and MDSeq ability to identify differentially dispersed genes.
figure_basename <- "fig1-DiPhiSeq_vs_MDSeq-DE-AUC"
## simulation parameters
### DE genes
DE_fraction <- 0.5

pdf(sprintf("%s/%s.pdf", output_dir, figure_basename), width=9)
for (DE_base_effect in c(1.1, 1.2, 1.3, 1.4, 1.5)) {
  DE_effect_param <- get_DE_effect_param(DE_base_effect)
  df2ggplot_all_genes_AUC <- data.frame()
  # read AUC statistics files
  for (one_sample_size in sample_sizes) {
    for (one_replicate in 1:replicates) {
      dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, one_sample_size, one_sample_size, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
      ## DiPhiSeq
      diphiseq_output_dir <- sprintf("%s/10-DiPhiSeq/00-highly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      diphiseq_output_basename <- sprintf("%s_%s_DiPhiSeq_results", dataset_name, normalization_method)
      auc_file <- sprintf("%s/%s_auc_stats.csv", diphiseq_output_dir, diphiseq_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("DiPhiSeq", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
      ## MDSeq
      mdseq_output_dir <- sprintf("%s/20-MDSeq/00-highly_DE/base_effect_%s/%s/repl_%d/00-all_genes", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      mdseq_output_basename <- sprintf("%s_%s_MDSeq_FC_%d", dataset_name, normalization_method, DD_FC_threshold)
      auc_file <- sprintf("%s/%s_auc_stats.csv", mdseq_output_dir, mdseq_output_basename)
      if (file.exists(auc_file)) {
        auc_df <- read.csv(file=auc_file)
        auc_df <- auc_df[which(auc_df$category %in% c("all", "DE", "nonDE")),]
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(auc_df)[1]), samp_cond_2=rep(one_sample_size, dim(auc_df)[1]), one_replicate=rep(one_replicate, dim(auc_df)[1]), method=rep("MDSeq", dim(auc_df)[1]), category=as.character(auc_df$category), value=auc_df$AUC)
        df2ggplot_all_genes_AUC <- rbind(df2ggplot_all_genes_AUC, one_df)
      }
    }
  }
  if (nrow(df2ggplot_all_genes_AUC) > 0) {
    df2ggplot_all_genes_AUC$samp_cond_1 <- as.factor(df2ggplot_all_genes_AUC$samp_cond_1)
    # plot
    plot_title <- sprintf("Differential dispersion performances - AUC - %d replicates \nDE: pct=%.1f, base effect=%.1f, effect param=%.2f, runif: %d;\nDV: pct=%.1f, base effect=%.1f, effect param=%.2f; runif: %d\nOutliers: single: high=%.2f, low=%.2f; random: high=%.2f, low=%.2f\nfilter: %d, outlier removal: %d\nMDSeq FC thresholds: mean=%.2f, dispersion=%.2f\nlabel FC threshold: mean=%.2f, dispersion=%.2f\np-val correction: DiPhiSeq:%s, MDSeq:%s", replicates, DE_fraction, DE_base_effect, DE_effect_param, non_DE_runif_FC, DD_fraction, DD_base_effect, DD_effect_param, non_DD_runif_FC, single_outlier_high_fraction, single_outlier_low_fraction, random_outlier_high_fraction, random_outlier_low_fraction, filter_threshold, outlier_removal_bool_param, DE_base_effect, DD_FC_threshold, DE_base_effect, DD_base_effect, diphiseq_pval_correction, mdseq_pval_correction)
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, plot_title, col = "black")
    p <- ggplot(data = df2ggplot_all_genes_AUC, aes(x=samp_cond_1, y=value, fill=category)) + 
      geom_boxplot() +
      facet_wrap(~method) +
      theme_bw() +
      theme(panel.border = element_rect(color="grey50")) +
      theme(legend.position = "bottom") +
      theme(legend.title=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16)) +
      theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=14), strip.text = element_text(size=14)) +
      labs(x = "Sample size per condition", y = "AUC", fill="Differential mean\nexpression category") +
      scale_fill_discrete(labels = c("all genes", "highly DE", "lowly DE"))
    print(p)
  }
}
dev.off()
write.csv(df2ggplot_all_genes_AUC, file=sprintf("%s/%s.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)


# Figure 2: DiPhiSeq and MDSeq ability to detect differential dispersion for lowly differentially expressed genes.
figure_basename <- "fig2-DiPhiSeq_vs_MDSeq-lowly_DE-FDR_TPR"
## simulation parameters
### DE genes
DE_fraction <- 0
DE_effect_param <- 0

df2ggplot <- data.frame()
for (DE_base_effect in c(1.1, 1.2, 1.3, 1.4, 1.5)) {
  for (one_sample_size in sample_sizes) {
    for (one_replicate in 1:replicates) {
      dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, one_sample_size, one_sample_size, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
      # read result performance files
      ## DiPhiSeq
      diphiseq_output_dir <- sprintf("%s/10-DiPhiSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      diphiseq_output_basename <- sprintf("%s_%s_DiPhiSeq_results", dataset_name, normalization_method)
      perf_file <- sprintf("%s/%s_performance_stats.csv", diphiseq_output_dir, diphiseq_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("DiPhiSeq", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
      ## MDSeq
      mdseq_output_dir <- sprintf("%s/20-MDSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      mdseq_output_basename <- sprintf("%s_%s_MDSeq_FC_%d", dataset_name, normalization_method, DD_FC_threshold)
      perf_file <- sprintf("%s/%s_performance_stats.csv", mdseq_output_dir, mdseq_output_basename)
      if (file.exists(perf_file)) {
        perf_df <- read.csv(file=perf_file)
        one_df <- data.frame(samp_cond_1=rep(one_sample_size, dim(perf_df)[1]), samp_cond_2=rep(one_sample_size, dim(perf_df)[1]), max_mean_FC=rep(DE_base_effect, dim(perf_df)[1]), replicate=rep(one_replicate, dim(perf_df)[1]), method=rep("MDSeq", dim(perf_df)[1]), stat=as.character(perf_df$stat), value=perf_df$value)
        df2ggplot <- rbind(df2ggplot, one_df)
      }
    }
  }
}

pdf(sprintf("%s/%s.pdf", output_dir, figure_basename), width=9)
# multiple plot
## FDR
one_stat <- "FDR"
one_df2ggplot <- df2ggplot[which(df2ggplot$stat==one_stat),]
one_df2ggplot$samp_cond_1 <- factor(one_df2ggplot$samp_cond_1)
one_df2ggplot$samp_cond_2 <- factor(one_df2ggplot$samp_cond_2)
one_df2ggplot$max_mean_FC <- factor(one_df2ggplot$max_mean_FC)
fdr_plot <- ggplot(data = one_df2ggplot, aes(x=samp_cond_1, y=value, fill=max_mean_FC)) + 
  geom_boxplot() +
  geom_hline(aes(yintercept=pval_threshold), linetype="dotted") +
  facet_wrap(~method) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(legend.title=element_text(size=14), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=12), strip.text = element_text(size=14)) +
  labs(x = "Sample size per condition", y = one_stat, fill="Maximum mean fold-change value")
### get the legend and remove it
plot_legend <- get_legend(fdr_plot)
fdr_plot <- fdr_plot + theme(legend.position = "none")
## TPR
one_stat <- "TPR"
one_df2ggplot <- df2ggplot[which(df2ggplot$stat==one_stat),]
one_df2ggplot$samp_cond_1 <- factor(one_df2ggplot$samp_cond_1)
one_df2ggplot$samp_cond_2 <- factor(one_df2ggplot$samp_cond_2)
one_df2ggplot$max_mean_FC <- factor(one_df2ggplot$max_mean_FC)
tpr_plot <- ggplot(data = one_df2ggplot, aes(x=samp_cond_1, y=value, fill=max_mean_FC)) + 
  geom_boxplot() +
  facet_wrap(~method) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), strip.text = element_text(size=14)) +
  labs(x = "Sample size per condition", y = one_stat)
## multiplot
multiplot_facet <- ggdraw() +
  draw_plot(fdr_plot, 0, 0.1, 0.5, 0.9) +
  draw_plot(tpr_plot, 0.5, 0.1, 0.5, 0.9) +
  draw_plot(plot_legend, 0, 0, 1, 0.1) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), size=15)
print(multiplot_facet)
dev.off()
write.csv(df2ggplot, file=sprintf("%s/%s.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)


# Figure 3: Truly DD genes identified by either DiPhiSeq, MDSeq, or both,among lowly differentially expressed genes.
figure_basename <- "fig3-DiPhiSeq_vs_MDSeq-truly_DD"
# path to output directory
output_dir <- "./output/test/simulations"
## simulation parameters
### DE genes
sample_sizes <- c(30, 40, 50, 100)
DE_fraction <- 0
DE_effect_param <- 0
DE_base_effects <- c(1.1, 1.2, 1.3, 1.4, 1.5)
diphiseq_DD_pval_colname <- "fdr.phi" # use p-values corrected by the Benjamini-Yekutieli FDR-controlling procedure
diphiseq_DD_pval_colname_category <- sprintf("category.%s", diphiseq_DD_pval_colname)
mdseq_DD_pval_colname <- "FDR.dispersion" # use p-values corrected by the Benjamini-Yekutieli FDR-controlling procedure
mdseq_DD_pval_colname_category <- sprintf("category.%s", mdseq_DD_pval_colname)

all_count_ggplot_df <- data.frame()
all_df2ggplot <- data.frame()
for (one_sample_size in sample_sizes) {
  print(sprintf("nb samples: %d", one_sample_size))
  for (DE_base_effect in DE_base_effects) {
    print(sprintf("base effect DE: %.1f", DE_base_effect))
    for (one_replicate in seq(1:replicates)) {
      # load simulated data to get true mean and dispersion fold-changes
      data_dir <- sprintf("%s/00-Data/10-lowly_DE/base_effect_%s/%d/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
      dataset_name <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s_repl_%d", genes, one_sample_size, one_sample_size, depth, gsub("[.]", "_", DE_fraction), gsub("[.]", "_", DE_base_effect), gsub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, gsub("[.]", "_", DD_fraction), gsub("[.]", "_", DD_base_effect), gsub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, gsub("[.]", "_", single_outlier_high_fraction), gsub("[.]", "_", single_outlier_low_fraction), gsub("[.]", "_", random_outlier_high_fraction), gsub("[.]", "_", random_outlier_low_fraction), one_replicate)
      rds_file <- sprintf("%s/%s.rds", data_dir, dataset_name)
      if (file.exists(rds_file)) {
        sim_data <- readRDS(rds_file) # contains compData object (compcodeR)
        ## simulation parameters
        annotations_df <- sim_data@variable.annotations
        ### add true log2FC for dispersion
        truelog2foldchanges.dispersion <- log(annotations_df$truedispersions.S2 / annotations_df$truedispersions.S1, 2)
        annotations_df <- cbind(annotations_df, truelog2foldchanges.dispersion)
        
        # read DiPhiSeq and MDSeq output tables
        ## DiPhiSeq
        diphiseq_output_dir <- sprintf("%s/10-DiPhiSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
        diphiseq_result_performance_file <- sprintf("%s/%s_%s_DiPhiSeq_results_performance.csv", diphiseq_output_dir, dataset_name, normalization_method)
        ## MDSeq
        mdseq_output_dir <- sprintf("%s/20-MDSeq/10-lowly_DE/base_effect_%s/%s/repl_%d", output_dir, sub("[.]", "_", DE_base_effect), one_sample_size, one_replicate)
        mdseq_result_performance_file <- sprintf("%s/%s_%s_MDSeq_FC_%d_performance.csv", mdseq_output_dir, dataset_name, normalization_method, DD_FC_threshold)
        if (file.exists(diphiseq_result_performance_file) & file.exists(mdseq_result_performance_file)) {
          diphiseq_result_performance <- read.csv(file=diphiseq_result_performance_file, row.names=1)
          mdseq_result_performance <- read.csv(file=mdseq_result_performance_file, row.names=1)
          print(sprintf("files for nb samples: %d and base effect DE: %.1f", one_sample_size, DE_base_effect))
          # get results per performance category
          for (one_category in c("TP", "FP", "TN", "FN")) {
            ## DiPhiSeq
            diphiseq_one_category_genes <- rownames(diphiseq_result_performance[which(diphiseq_result_performance[[diphiseq_DD_pval_colname_category]] == one_category),])
            diphiseq_one_category_gene_count <- length(diphiseq_one_category_genes)
            ## MDSeq
            mdseq_one_category_genes <- rownames(mdseq_result_performance[which(mdseq_result_performance[[mdseq_DD_pval_colname_category]] == one_category),])
            mdseq_one_category_gene_count <- length(mdseq_one_category_genes)
            ## common genes
            one_category_common_genes <- diphiseq_one_category_genes[diphiseq_one_category_genes %in% mdseq_one_category_genes]
            one_category_common_gene_count <- length(one_category_common_genes)
            ## DiPhiSeq specific genes
            diphiseq_one_category_specific_genes <- diphiseq_one_category_genes[! diphiseq_one_category_genes %in% one_category_common_genes]
            diphiseq_one_category_specific_gene_count <- length(diphiseq_one_category_specific_genes)
            ## MDSeq specific genes
            mdseq_one_category_specific_genes <- mdseq_one_category_genes[! mdseq_one_category_genes %in% one_category_common_genes]
            mdseq_one_category_specific_gene_count <- length(mdseq_one_category_specific_genes)
            ## total category genes
            one_category_total_gene_count <- one_category_common_gene_count + diphiseq_one_category_specific_gene_count + mdseq_one_category_specific_gene_count
            
            # all count ggplot data frame
            one_category_count_vector <- c(diphiseq_one_category_gene_count, one_category_common_gene_count, diphiseq_one_category_specific_gene_count, mdseq_one_category_gene_count, one_category_common_gene_count, mdseq_one_category_specific_gene_count)
            all_one_category_counts <- length(one_category_count_vector)
            category_vector <- rep(one_category, all_one_category_counts)
            method_vector <- c(rep("DiPhiSeq", all_one_category_counts/2), rep("MDSeq", all_one_category_counts/2))
            venn_category_vector <- rep(c("all", "common", "specific"), 2)
            count_ggplot_df <- data.frame(sample_size_cond_1=rep(one_sample_size, all_one_category_counts), 
                                          sample_size_cond_2=rep(one_sample_size, all_one_category_counts), 
                                          DE_base_effect=rep(DE_base_effect, all_one_category_counts), 
                                          replicate=rep(one_replicate, all_one_category_counts), 
                                          category=category_vector, 
                                          method=method_vector, 
                                          venn_category=venn_category_vector, 
                                          count=one_category_count_vector,
                                          pct=one_category_count_vector/one_category_total_gene_count*100)
            all_count_ggplot_df <- rbind(all_count_ggplot_df, count_ggplot_df)
            
            # density ggplot data frame
            annotations_sub_df <- annotations_df[rownames(annotations_df) %in% diphiseq_one_category_specific_genes | rownames(annotations_df) %in% mdseq_one_category_specific_genes | rownames(annotations_df) %in% one_category_common_genes, ]
            nb_genes <- length(rownames(annotations_sub_df))
            ## venn categories
            venn_category_vector_2 <- unlist(lapply(rownames(annotations_sub_df), common=one_category_common_genes, diphiseq=diphiseq_one_category_specific_genes, mdseq=mdseq_one_category_specific_genes, function(x, common, diphiseq, mdseq) {
              if (x %in% common) {
                return("DiPhiSeq and MDSeq")
              } else {
                if (x %in% diphiseq) {
                  return("DiPhiSeq only")
                } else {
                  if (x %in% mdseq) {
                    return("MDSeq only")
                  }
                }
              }
            }))
            mean_exp <- log2(unlist(lapply(rownames(annotations_sub_df), df=annotations_sub_df, function(x, df) { return(mean(df[x, "truemeans.S1"] + df[x, "truemeans.S2"]))})))
            
            df2ggplot <- data.frame(sample_size_cond_1=rep(one_sample_size, nb_genes), 
                                    sample_size_cond_2=rep(one_sample_size, nb_genes), 
                                    DE_base_effect=rep(DE_base_effect, nb_genes), 
                                    replicate=rep(one_replicate, nb_genes),
                                    gene=rownames(annotations_sub_df), 
                                    category=rep(one_category, nb_genes), 
                                    venn.category=venn_category_vector_2, 
                                    true_FC_mean=annotations_sub_df[rownames(annotations_sub_df), "truelog2foldchanges"], 
                                    true_FC_dispersion=annotations_sub_df[rownames(annotations_sub_df), "truelog2foldchanges.dispersion"], 
                                    mean=mean_exp)
            all_df2ggplot <- rbind(all_df2ggplot, df2ggplot)
          }
        }
      }
    }
  }
}
# write csv files
write.csv(all_count_ggplot_df, file=sprintf("%s/%s_category_counts.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)
write.csv(all_df2ggplot, file=sprintf("%s/%s_category_true_FCs_mean_exp.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)

# compute mean counts for TP
all_TP_mean_count_ggplot_df <- data.frame()
TP_count_ggplot_df <- all_count_ggplot_df[which(all_count_ggplot_df$category == "TP"),]
methods <- c("DiPhiSeq", "MDSeq")
venn_categories <- c("all", "common", "specific")
sample_sizes <- unique(TP_count_ggplot_df$sample_size_cond_1)
DE_base_effects <- unique(TP_count_ggplot_df$DE_base_effect)
for (one_sample_size in sample_sizes) {
  print(sprintf("nb samples: %d", one_sample_size))
  for (one_DE_base_effect in DE_base_effects) {
    print(sprintf("base effect DE: %.1f", DE_base_effect))
    sample_size_vector <- DE_base_effect_vector <- method_vector <- venn_category_vector <- venn_category_vector <- mean_count_vector <- c()
    for (one_method in methods) {
      for (one_venn_category in venn_categories) {
        mean_count <- round(mean(TP_count_ggplot_df[which(TP_count_ggplot_df$sample_size_cond_1==one_sample_size & TP_count_ggplot_df$sample_size_cond_2==one_sample_size & TP_count_ggplot_df$DE_base_effect==one_DE_base_effect & TP_count_ggplot_df$method==one_method & TP_count_ggplot_df$venn_category==one_venn_category), "count"]))
        mean_count_vector <- c(mean_count_vector, mean_count)
        venn_category_vector <- c(venn_category_vector, one_venn_category)
      }
      method_vector <- c(method_vector, rep(one_method, length(venn_categories)))
    }
    
    TP_mean_count_ggplot_df <- data.frame(method=method_vector,
                                          venn_category=venn_category_vector,
                                          mean_count=mean_count_vector)
    ## add percentages
    total_gene_count <- TP_mean_count_ggplot_df[which(TP_mean_count_ggplot_df$method=="DiPhiSeq" & TP_mean_count_ggplot_df$venn_category=="common"), "mean_count"] + TP_mean_count_ggplot_df[which(TP_mean_count_ggplot_df$method=="DiPhiSeq" & TP_mean_count_ggplot_df$venn_category=="specific"), "mean_count"] + TP_mean_count_ggplot_df[which(TP_mean_count_ggplot_df$method=="MDSeq" & TP_mean_count_ggplot_df$venn_category=="specific"), "mean_count"]
    pcts <- TP_mean_count_ggplot_df[, "mean_count"] / total_gene_count * 100
    TP_mean_count_ggplot_df <- cbind(TP_mean_count_ggplot_df, pct=pcts)
    ## add sample size and DE base effect
    reps <- dim(TP_mean_count_ggplot_df)[1]
    TP_mean_count_ggplot_df <- cbind(sample_size_cond_1=rep(one_sample_size, reps), samp_size_cond_2=rep(one_sample_size, reps), DE_base_effect=rep(one_DE_base_effect, reps), TP_mean_count_ggplot_df)
    all_TP_mean_count_ggplot_df <- rbind(all_TP_mean_count_ggplot_df, TP_mean_count_ggplot_df)
  }
}
# write csv files
write.csv(all_TP_mean_count_ggplot_df, file=sprintf("%s/%s_TP_mean_counts.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)

# plots
sample_size_to_plot <- 50
DE_base_effect_to_plot <- 1.5
replicates <- 3

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # no log file
## mean Venn diagrams
one_sample_size_one_DE_base_effect_TP_mean_count_df <- all_TP_mean_count_ggplot_df[which(all_TP_mean_count_ggplot_df$sample_size_cond_1==sample_size_to_plot & all_TP_mean_count_ggplot_df$samp_size_cond_2==sample_size_to_plot & all_TP_mean_count_ggplot_df$DE_base_effect==DE_base_effect_to_plot),]
diphiseq_circle_area <- one_sample_size_one_DE_base_effect_TP_mean_count_df[which(one_sample_size_one_DE_base_effect_TP_mean_count_df$method == "DiPhiSeq" & one_sample_size_one_DE_base_effect_TP_mean_count_df$venn_category == "all"), "mean_count"]
mdseq_circle_area <- one_sample_size_one_DE_base_effect_TP_mean_count_df[which(one_sample_size_one_DE_base_effect_TP_mean_count_df$method == "MDSeq" & one_sample_size_one_DE_base_effect_TP_mean_count_df$venn_category == "all"), "mean_count"]
common_area <- one_sample_size_one_DE_base_effect_TP_mean_count_df[which(one_sample_size_one_DE_base_effect_TP_mean_count_df$method == "MDSeq" & one_sample_size_one_DE_base_effect_TP_mean_count_df$venn_category == "common"), "mean_count"]
venn_plot <- venn.diagram(x=NULL, filename=NULL, direct.area=TRUE, area.vector=c(diphiseq_circle_area, mdseq_circle_area, common_area), category=c("DiPhiSeq", "MDSeq"), fill=c(brewer.pal(9, "Set1")[1], brewer.pal(9, "Set1")[2]), print.mode=c("raw", "percent"), sigdigs=3) 

## mean expression distributions
one_sample_size_one_DE_base_effect_TP_df2ggplot <- all_df2ggplot[which(all_df2ggplot$sample_size_cond_1 == sample_size_to_plot & all_df2ggplot$sample_size_cond_2 == sample_size_to_plot & all_df2ggplot$DE_base_effect == DE_base_effect_to_plot & all_df2ggplot$category == "TP"),]
### density facets
df2ggplot_density_facet <- data.frame(
  value=c(abs(one_sample_size_one_DE_base_effect_TP_df2ggplot[, "true_FC_mean"]), abs(one_sample_size_one_DE_base_effect_TP_df2ggplot[, "true_FC_dispersion"]), one_sample_size_one_DE_base_effect_TP_df2ggplot[, "mean"]),
  stat=c(rep("True mean FC", length(one_sample_size_one_DE_base_effect_TP_df2ggplot[, "true_FC_mean"])), rep("True dispersion FC", length(one_sample_size_one_DE_base_effect_TP_df2ggplot[, "true_FC_dispersion"])), rep("Mean expression", length(one_sample_size_one_DE_base_effect_TP_df2ggplot[, "mean"]))),
  venn.category=rep(one_sample_size_one_DE_base_effect_TP_df2ggplot[, "venn.category"], 3)
)
df2ggplot_density_facet$venn.category <- factor(df2ggplot_density_facet$venn.category, levels=c("DiPhiSeq and MDSeq", "DiPhiSeq only", "MDSeq only"))
levels(df2ggplot_density_facet$venn.category)[which(levels(df2ggplot_density_facet$venn.category) == "DiPhiSeq and MDSeq")] <- "DiPhiSeq\nand MDSeq"
### t-tests and facet titles
values <- levels(df2ggplot_density_facet$stat)
facet_titles <- unlist(lapply(values, df=df2ggplot_density_facet, function(x, df) {
  diphiseq_specific_values <- df[which(df$stat==x & df$venn.category=="DiPhiSeq only"), "value"]
  mdseq_specific_values <- df[which(df$stat==x & df$venn.category=="MDSeq only"), "value"]
  if (x == "True dispersion FC") {
    test_result <- wilcox.test(diphiseq_specific_values, mdseq_specific_values, alternative="greater")
  } else {
    test_result <- wilcox.test(mdseq_specific_values, diphiseq_specific_values, alternative="greater")
  }
  return(sprintf("%s\np-value = %.2e", x, test_result$p.value))
}))
names(facet_titles) = values

density_facet.plot <- ggplot(df2ggplot_density_facet, aes(x=value, col=venn.category)) +
  geom_line(stat="density") +
  facet_wrap(~stat, scales="free", labeller=labeller(stat=facet_titles)) +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Value", y="Density") +
  theme(legend.title=element_blank()) +
  scale_color_manual(values=c(brewer.pal(9, "Set1")[4], brewer.pal(9, "Set1")[1], brewer.pal(9, "Set1")[2]))

# venn categories barplot for all sample sizes
venn_barplot_ggplot_df <- data.frame()
sample_sizes <- unique(all_TP_mean_count_ggplot_df$sample_size_cond_1)
DE_base_effects <- unique(all_TP_mean_count_ggplot_df$DE_base_effect)
for (one_sample_size in sample_sizes) {
  print(sprintf("nb samples: %d", one_sample_size))
  for (one_DE_base_effect in DE_base_effects) {
    print(sprintf("base effect DE: %.1f", DE_base_effect))
    sample_size_vector <- DE_base_effect_vector <- method_vector <- venn_category_vector <- venn_category_vector <- mean_count_vector <- c()
    one_sample_size_one_DE_base_effect_TP_mean_count_df <- all_TP_mean_count_ggplot_df[which(all_TP_mean_count_ggplot_df$sample_size_cond_1==one_sample_size & all_TP_mean_count_ggplot_df$samp_size_cond_2==one_sample_size & all_TP_mean_count_ggplot_df$DE_base_effect==one_DE_base_effect),]
    ## common
    common_mean_count_df <- one_sample_size_one_DE_base_effect_TP_mean_count_df[which(one_sample_size_one_DE_base_effect_TP_mean_count_df$method=="DiPhiSeq" & one_sample_size_one_DE_base_effect_TP_mean_count_df$venn_category=="common"),]
    common_mean_count <- common_mean_count_df[, "mean_count"]
    common_pct <- common_mean_count_df[, "pct"]
    ## DiPhiSeq specific
    diphiseq_specific_mean_count_df <- one_sample_size_one_DE_base_effect_TP_mean_count_df[which(one_sample_size_one_DE_base_effect_TP_mean_count_df$method=="DiPhiSeq" & one_sample_size_one_DE_base_effect_TP_mean_count_df$venn_category=="specific"),]
    diphiseq_specific_mean_count <- diphiseq_specific_mean_count_df[, "mean_count"]
    diphiseq_specific_pct <- diphiseq_specific_mean_count_df[, "pct"]
    ## MDSeq specific
    mdseq_specific_mean_count_df <- one_sample_size_one_DE_base_effect_TP_mean_count_df[which(one_sample_size_one_DE_base_effect_TP_mean_count_df$method=="MDSeq" & one_sample_size_one_DE_base_effect_TP_mean_count_df$venn_category=="specific"),]
    mdseq_specific_mean_count <- mdseq_specific_mean_count_df[, "mean_count"]
    mdseq_specific_pct <- mdseq_specific_mean_count_df[, "pct"]
    
    one_sample_size_one_DE_base_effect_TP_venn_df <- cbind(sample_size_cond_1=rep(one_sample_size, 3), samp_size_cond_2=rep(one_sample_size, 3), DE_base_effect=rep(one_DE_base_effect, 3), category=c("DiPhiSeq\nand MDSeq", "DiPhiSeq only", "MDSeq only"), mean_count=c(common_mean_count, diphiseq_specific_mean_count, mdseq_specific_mean_count), pct=c(common_pct, diphiseq_specific_pct, mdseq_specific_pct))
    venn_barplot_ggplot_df <- rbind(venn_barplot_ggplot_df, one_sample_size_one_DE_base_effect_TP_venn_df)
  }
}
venn_barplot_ggplot_df$mean_count <- as.numeric(as.character(venn_barplot_ggplot_df$mean_count))
venn_barplot_ggplot_df$pct <- as.numeric(as.character(venn_barplot_ggplot_df$pct))
venn_barplot_ggplot_df$sample_size_cond_1 <- factor(venn_barplot_ggplot_df$sample_size_cond_1, levels=rev(levels(venn_barplot_ggplot_df$sample_size_cond_1)))
venn_barplot_ggplot_df$category <- factor(venn_barplot_ggplot_df$category, levels=c("MDSeq only", "DiPhiSeq\nand MDSeq", "DiPhiSeq only"))
write.csv(venn_barplot_ggplot_df, file=sprintf("%s/%s_TP_mean_count_Venn_barplot.csv", output_dir, figure_basename), quote=FALSE, row.names=FALSE)

venn_barplot <- ggplot(venn_barplot_ggplot_df[which(venn_barplot_ggplot_df$DE_base_effect==DE_base_effect_to_plot),], aes(x=sample_size_cond_1, y=pct, fill=category)) +
  geom_bar(stat="identity") +
  coord_flip() +
  geom_text(aes(label=mean_count), size=4, position=position_stack(vjust=0.5)) +
  scale_fill_manual(values=c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[4], brewer.pal(9, "Set1")[1])) +
  labs(x="Sample size", y="Percentage") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_blank()) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10)) +
  theme(axis.ticks.y=element_blank()) +
  guides(fill=guide_legend(reverse=TRUE))

pdf(sprintf("%s/%s.pdf", output_dir, figure_basename), width=9)
multiplot_facet <- ggdraw() +
  draw_plot(venn_plot, 0.05, 0.625, 0.25, 0.25) +
  draw_plot(density_facet.plot, 0.35, 0.5, 0.65, 0.5) +
  draw_plot(venn_barplot, 0, 0, 1, 0.5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size=15)
print(multiplot_facet)
dev.off()

