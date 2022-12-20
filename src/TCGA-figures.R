library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(UpSetR)
library(grid)


##############
# Parameters #
##############
# datasets
datasets <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-THCA")
# normalization method
normalization_method <- "TMM"
# minimum expression threshold in CPM to filter out lowly expressed genes
filter_threshold <- 1
# MDSeq parameters
## perform MDSeq outlier removal function
outlier_removal <- TRUE
## fold-change threshold to identify DE genes (not in log2 scale)
DE_FC_threshold <- 1
## fold-change threshold to identify DD genes (not in log2 scale)
DD_FC_threshold <- 1
## p-value threshold
pval_threshold <- 0.05
# path to output directory
output_dir <- "./output/TCGA"
# first condition
pheno_1 <- "TP"
# second condition
pheno_2 <- "NT"
# method order
method_order <- c("Levene", "MDSeq", "DiPhiSeq", "GAMLSS", "DiffDist")


############
# Figure 4 #
############
fig_output_basename <- sprintf("fig4-TCGA_%s_vs_%s_FC_%d_DE_DD_counts", pheno_1, pheno_2, DD_FC_threshold)

# MDSeq results
mdseq_out_dir <- sprintf("%s/20-MDSeq", output_dir)
DE_vector <- DEp_vector <- DEm_vector <- non_DE_vector <- DVp_non_DE_vector <- DVm_non_DE_vector <- non_DE_non_DV_vector <- total_vector <- c()
all_datasets_MDSeq_genes <- list()
for (one_dataset in datasets) {
  dataset_mdseq_out_dir <- sprintf("%s/%s", mdseq_out_dir, one_dataset)
  ## all gene analysis
  ### read MDSeq output table
  mdseq_output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq_FC_%s", one_dataset, pheno_1, pheno_2, normalization_method, sub("[.]", "_", log(DD_FC_threshold, 2)))
  mdseq_results_file <- sprintf("%s/%s/%s.csv", mdseq_out_dir, one_dataset, mdseq_output_basename)
  
  if (file.exists(mdseq_results_file)) {
    mdseq_results_df <- read.csv(file=mdseq_results_file, header=TRUE, row.names=1)
    
    # total
    nb_total <- dim(mdseq_results_df)[1]
    
    # DE
    DE <- sum(mdseq_results_df$FDR.mean < pval_threshold, na.rm=TRUE)
    ## DE+
    DEp <- sum(mdseq_results_df$FDR.mean < pval_threshold & mdseq_results_df[[sprintf("%svs%s.mean.log2FC.%s", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ## DE-
    DEm <- sum(mdseq_results_df$FDR.mean < pval_threshold & mdseq_results_df[[sprintf("%svs%s.mean.log2FC.%s", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] < 0, na.rm=TRUE)
    # non-DE
    non_DE <- sum(mdseq_results_df$FDR.mean > pval_threshold, na.rm=TRUE)
    
    # get non-DE genes
    non_DE_genes <- rownames(mdseq_results_df[which(mdseq_results_df$FDR.mean > pval_threshold),])
    all_datasets_MDSeq_genes[[one_dataset]] <- list(non_DE=non_DE_genes)
    
    
    mdseq_dv_df <- mdseq_results_df[which(mdseq_results_df$FDR.mean > pval_threshold),]
    
    ## DV+ non-DE
    DVp_non_DE <- sum(mdseq_dv_df$FDR.dispersion < pval_threshold & mdseq_dv_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] > 0, na.rm=TRUE)
    DVp_non_DE_genes <- rownames(mdseq_dv_df[which(mdseq_dv_df$FDR.dispersion < pval_threshold & mdseq_dv_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] > 0),])
    all_datasets_MDSeq_genes[[one_dataset]][["DVp_non_DE"]] <- DVp_non_DE_genes
    
    ## DV- non-DE
    DVm_non_DE <- sum(mdseq_dv_df$FDR.dispersion < pval_threshold & mdseq_dv_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] < 0, na.rm=TRUE)
    DVm_non_DE_genes <- rownames(mdseq_dv_df[which(mdseq_dv_df$FDR.dispersion < pval_threshold & mdseq_dv_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] < 0),])
    all_datasets_MDSeq_genes[[one_dataset]][["DVm_non_DE"]] <- DVm_non_DE_genes
    
    ## non-DV non-DE
    non_DE_no_DV <- sum(mdseq_dv_df$FDR.dispersion > pval_threshold, na.rm=TRUE)
  } else {
    DE <- DEp <- DEm <- non_DE <- nb_total <-NA
    DVp_non_DE <- DVm_non_DE <- non_DE_no_DV <- NA
  }
  # concatenate vectors
  DE_vector <- c(DE_vector, DE)
  DEp_vector <- c(DEp_vector, DEp)
  DEm_vector <- c(DEm_vector, DEm)
  non_DE_vector <- c(non_DE_vector, non_DE)
  total_vector <- c(total_vector, nb_total)
  DVp_non_DE_vector <- c(DVp_non_DE_vector, DVp_non_DE)
  DVm_non_DE_vector <- c(DVm_non_DE_vector, DVm_non_DE)
  non_DE_non_DV_vector <- c(non_DE_non_DV_vector, non_DE_no_DV)
}
all_datasets_mdseq_df <- rbind(DE_vector, DEp_vector, DEm_vector, non_DE_vector, DVp_non_DE_vector, DVm_non_DE_vector, non_DE_non_DV_vector, total_vector)
colnames(all_datasets_mdseq_df) <- datasets
rownames(all_datasets_mdseq_df) <- c("DE", "DE+", "DE-", "non-DE", "DV+ non-DE", "DV- non-DE", "non-DV non-DE", "total")
write.csv(all_datasets_mdseq_df, file=sprintf("%s/%s_MDSeq.csv", output_dir, fig_output_basename), quote=FALSE, row.names=TRUE) 

# Levene's test results
Levene_out_dir <- sprintf("%s/10-Levene", output_dir)
all_datasets_Levene_DV_MDSeq_no_DE_genes <- list()
DVp_MDSeq_no_DE_vector <- DVm_MDSeq_no_DE_vector <- no_DV_MDSeq_no_DE_vector <- c()
for (one_dataset in datasets) {
  dataset_Levene_out_dir <- sprintf("%s/%s", Levene_out_dir, one_dataset)
  Levene_output_basename <- sprintf("%s_%s_vs_%s_%s_Levene_test_results", one_dataset, pheno_1, pheno_2, normalization_method)
  Levene_results_file <- sprintf("%s/%s.csv", dataset_Levene_out_dir, Levene_output_basename)
  
  if (file.exists(Levene_results_file)) {
    Levene_results_df <- read.csv(file=Levene_results_file, header=TRUE, row.names=1)
    
    ## DV+
    DVp <- sum(Levene_results_df$p_value.BH < pval_threshold & Levene_results_df[, sprintf("var_%s", pheno_1)] > Levene_results_df[, sprintf("var_%s", pheno_2)], na.rm=TRUE)
    DVp_genes <- rownames(Levene_results_df[which(Levene_results_df$p_value.BH < pval_threshold & Levene_results_df[, sprintf("var_%s", pheno_1)] > Levene_results_df[, sprintf("var_%s", pheno_2)]),])
    DVp_MDSeq_no_DE_genes <- DVp_genes[DVp_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DVp_MDSeq_no_DE <- length(DVp_MDSeq_no_DE_genes)
    all_datasets_Levene_DV_MDSeq_no_DE_genes[[one_dataset]][["DVp_MDSeq_no_DE"]] <- DVp_MDSeq_no_DE_genes
    
    ## DV-
    DVm <- sum(Levene_results_df$p_value.BH < pval_threshold & Levene_results_df[, sprintf("var_%s", pheno_1)] < Levene_results_df[, sprintf("var_%s", pheno_2)], na.rm=TRUE)
    DVm_genes <- rownames(Levene_results_df[which(Levene_results_df$p_value.BH < pval_threshold & Levene_results_df[, sprintf("var_%s", pheno_1)] < Levene_results_df[, sprintf("var_%s", pheno_2)]),])
    DVm_MDSeq_no_DE_genes <- DVm_genes[DVm_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DVm_MDSeq_no_DE <- length(DVm_MDSeq_no_DE_genes)
    all_datasets_Levene_DV_MDSeq_no_DE_genes[[one_dataset]][["DVm_MDSeq_no_DE"]] <- DVm_MDSeq_no_DE
    
    ## non-DV
    no_DV <- sum(Levene_results_df$p_value.BH > pval_threshold, na.rm=TRUE)
    no_DV_genes <- rownames(Levene_results_df[which(Levene_results_df$p_value.BH > pval_threshold),])
    no_DV_MDSeq_no_DE_genes <- no_DV_genes[no_DV_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    no_DV_MDSeq_no_DE <- length(no_DV_MDSeq_no_DE_genes)
    all_datasets_Levene_DV_MDSeq_no_DE_genes[[one_dataset]][["no_DV_MDSeq_no_DE"]] <- no_DV_MDSeq_no_DE
    
    
  } else {
    DVp_MDSeq_no_DE <- DVm_MDSeq_no_DE <- no_DV_MDSeq_no_DE <- NA
  }
  # concatenate vectors
  DVp_MDSeq_no_DE_vector <- c(DVp_MDSeq_no_DE_vector, DVp_MDSeq_no_DE)
  DVm_MDSeq_no_DE_vector <- c(DVm_MDSeq_no_DE_vector, DVm_MDSeq_no_DE)
  no_DV_MDSeq_no_DE_vector <- c(no_DV_MDSeq_no_DE_vector, no_DV_MDSeq_no_DE)
}
all_datasets_Levene_df <- rbind(DVp_MDSeq_no_DE_vector, DVm_MDSeq_no_DE_vector, no_DV_MDSeq_no_DE_vector)
colnames(all_datasets_Levene_df) <- datasets
rownames(all_datasets_Levene_df) <- c("DV+ MDSeq non-DE", "DV- MDSeq non-DE", "non-DV MDSeq non-DE")
write.csv(all_datasets_Levene_df, file=sprintf("%s/%s_Levene_test_DV_MDSeq_non_DE.csv", output_dir, fig_output_basename), quote=FALSE, row.names=TRUE) 


# DiPhiSeq results
diphiseq_out_dir <- sprintf("%s/30-DiPhiSeq", output_dir)
all_datasets_diphiseq_DD_MDSeq_non_DE_genes <- list()
DDp_MDSeq_no_DE_vector <- DDm_MDSeq_no_DE_vector <- no_DD_MDSeq_no_DE_vector <- c()
for (one_dataset in datasets) {
  dataset_diphiseq_out_dir <- sprintf("%s/%s", diphiseq_out_dir, one_dataset)
  diphiseq_output_basename <- sprintf("%s_%s_vs_%s_%s_DiPhiSeq_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diphiseq_results_file <- sprintf("%s/%s.csv", dataset_diphiseq_out_dir, diphiseq_output_basename)
  
  if (file.exists(diphiseq_results_file)) {
    diphiseq_results_df <- read.csv(file=diphiseq_results_file, header=TRUE, row.names=1)
    
    ## DD+
    DDp <- sum(diphiseq_results_df$fdr.phi < pval_threshold & diphiseq_results_df[, sprintf("phi.%s", pheno_1)] > diphiseq_results_df[, sprintf("phi.%s", pheno_2)], na.rm=TRUE)
    DDp_genes <- rownames(diphiseq_results_df[which(diphiseq_results_df$fdr.phi < pval_threshold & diphiseq_results_df[, sprintf("phi.%s", pheno_1)] > diphiseq_results_df[, sprintf("phi.%s", pheno_2)]),])
    DDp_MDSeq_no_DE_genes <- DDp_genes[DDp_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DDp_MDSeq_no_DE <- length(DDp_MDSeq_no_DE_genes)
    all_datasets_diphiseq_DD_MDSeq_non_DE_genes[[one_dataset]][["DDp_MDSeq_no_DE"]] <- DDp_MDSeq_no_DE_genes
    
    ## DD-
    DDm <- sum(diphiseq_results_df$fdr.phi < pval_threshold & diphiseq_results_df[, sprintf("phi.%s", pheno_1)] < diphiseq_results_df[, sprintf("phi.%s", pheno_2)], na.rm=TRUE)
    DDm_genes <- rownames(diphiseq_results_df[which(diphiseq_results_df$fdr.phi < pval_threshold & diphiseq_results_df[, sprintf("phi.%s", pheno_1)] < diphiseq_results_df[, sprintf("phi.%s", pheno_2)]),])
    DDm_MDSeq_no_DE_genes <- DDm_genes[DDm_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DDm_MDSeq_no_DE <- length(DDm_MDSeq_no_DE_genes)
    all_datasets_diphiseq_DD_MDSeq_non_DE_genes[[one_dataset]][["DDm_MDSeq_no_DE"]] <- DDm_MDSeq_no_DE
    
    ## non-DD
    no_DD <- sum(diphiseq_results_df$fdr.phi > pval_threshold, na.rm=TRUE)
    no_DD_genes <- rownames(diphiseq_results_df[which(diphiseq_results_df$fdr.phi > pval_threshold),])
    no_DD_MDSeq_no_DE_genes <- no_DD_genes[no_DD_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    no_DD_MDSeq_no_DE <- length(no_DD_MDSeq_no_DE_genes)
    all_datasets_diphiseq_DD_MDSeq_non_DE_genes[[one_dataset]][["no_DD_MDSeq_no_DE"]] <- no_DD_MDSeq_no_DE
  } else {
    DDp_MDSeq_no_DE <- DDm_MDSeq_no_DE <- no_DD_MDSeq_no_DE <- NA
  }
  # concatenate vectors
  DDp_MDSeq_no_DE_vector <- c(DDp_MDSeq_no_DE_vector, DDp_MDSeq_no_DE)
  DDm_MDSeq_no_DE_vector <- c(DDm_MDSeq_no_DE_vector, DDm_MDSeq_no_DE)
  no_DD_MDSeq_no_DE_vector <- c(no_DD_MDSeq_no_DE_vector, no_DD_MDSeq_no_DE)
}
all_datasets_diphiseq_df <- rbind(DDp_MDSeq_no_DE_vector, DDm_MDSeq_no_DE_vector, no_DD_MDSeq_no_DE_vector)
colnames(all_datasets_diphiseq_df) <- datasets
rownames(all_datasets_diphiseq_df) <- c("DD+ MDSeq non-DE", "DD- MDSeq non-DE", "non-DD MDSeq non-DE")
write.csv(all_datasets_diphiseq_df, file=sprintf("%s/%s_DiPhiSeq_DD_MDSeq_non_DE.csv", output_dir, fig_output_basename), quote=FALSE, row.names=TRUE) 

# GAMLSS results
gamlss_out_dir <- sprintf("%s/40-GAMLSS", output_dir)
all_datasets_gamlss_DD_MDSeq_non_DE_genes <- list()
DDp_MDSeq_no_DE_vector <- DDm_MDSeq_no_DE_vector <- no_DD_MDSeq_no_DE_vector <- c()
for (one_dataset in datasets) {
  dataset_gamlss_out_dir <- sprintf("%s/%s", gamlss_out_dir, one_dataset)
  gamlss_output_basename <- sprintf("%s_%s_vs_%s_%s_GAMLSS_results", one_dataset, pheno_1, pheno_2, normalization_method)
  gamlss_results_file <- sprintf("%s/%s.csv", dataset_gamlss_out_dir, gamlss_output_basename)
  
  if (file.exists(gamlss_results_file)) {
    gamlss_results_df <- read.csv(file=gamlss_results_file, header=TRUE, row.names=1)
    
    ## DD+
    DDp <- sum(gamlss_results_df$padj.cv < pval_threshold & gamlss_results_df[, sprintf("CV.%s", pheno_1)] > gamlss_results_df[, sprintf("CV.%s", pheno_2)], na.rm=TRUE)
    DDp_genes <- rownames(gamlss_results_df[which(gamlss_results_df$padj.cv < pval_threshold & gamlss_results_df[, sprintf("CV.%s", pheno_1)] > gamlss_results_df[, sprintf("CV.%s", pheno_2)]),])
    DDp_MDSeq_no_DE_genes <- DDp_genes[DDp_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DDp_MDSeq_no_DE <- length(DDp_MDSeq_no_DE_genes)
    all_datasets_gamlss_DD_MDSeq_non_DE_genes[[one_dataset]][["DDp_MDSeq_no_DE"]] <- DDp_MDSeq_no_DE_genes
    
    ## DD-
    DDm <- sum(gamlss_results_df$padj.cv < pval_threshold & gamlss_results_df[, sprintf("CV.%s", pheno_1)] < gamlss_results_df[, sprintf("CV.%s", pheno_2)], na.rm=TRUE)
    DDm_genes <- rownames(gamlss_results_df[which(gamlss_results_df$padj.cv < pval_threshold & gamlss_results_df[, sprintf("CV.%s", pheno_1)] < gamlss_results_df[, sprintf("CV.%s", pheno_2)]),])
    DDm_MDSeq_no_DE_genes <- DDm_genes[DDm_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DDm_MDSeq_no_DE <- length(DDm_MDSeq_no_DE_genes)
    all_datasets_gamlss_DD_MDSeq_non_DE_genes[[one_dataset]][["DDm_MDSeq_no_DE"]] <- DDm_MDSeq_no_DE
    
    ## non-DD
    no_DD <- sum(gamlss_results_df$padj.cv > pval_threshold, na.rm=TRUE)
    no_DD_genes <- rownames(gamlss_results_df[which(gamlss_results_df$padj.cv > pval_threshold),])
    no_DD_MDSeq_no_DE_genes <- no_DD_genes[no_DD_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    no_DD_MDSeq_no_DE <- length(no_DD_MDSeq_no_DE_genes)
    all_datasets_gamlss_DD_MDSeq_non_DE_genes[[one_dataset]][["no_DD_MDSeq_no_DE"]] <- no_DD_MDSeq_no_DE
    
    
  } else {
    DDp_MDSeq_no_DE <- DDm_MDSeq_no_DE <- no_DD_MDSeq_no_DE <- NA
  }
  # concatenate vectors
  DDp_MDSeq_no_DE_vector <- c(DDp_MDSeq_no_DE_vector, DDp_MDSeq_no_DE)
  DDm_MDSeq_no_DE_vector <- c(DDm_MDSeq_no_DE_vector, DDm_MDSeq_no_DE)
  no_DD_MDSeq_no_DE_vector <- c(no_DD_MDSeq_no_DE_vector, no_DD_MDSeq_no_DE)
}
all_datasets_gamlss_df <- rbind(DDp_MDSeq_no_DE_vector, DDm_MDSeq_no_DE_vector, no_DD_MDSeq_no_DE_vector)
colnames(all_datasets_gamlss_df) <- datasets
rownames(all_datasets_gamlss_df) <- c("DD+ MDSeq non-DE", "DD- MDSeq non-DE", "non-DD MDSeq non-DE")
write.csv(all_datasets_gamlss_df, file=sprintf("%s/%s_GAMLSS_DD_MDSeq_non_DE.csv", output_dir, fig_output_basename), quote=FALSE, row.names=TRUE) 

# DiffDist results
diffdist_out_dir <- sprintf("%s/50-DiffDist", output_dir)
all_datasets_diffdist_DD_MDSeq_non_DE_genes <- list()
DDp_MDSeq_no_DE_vector <- DDm_MDSeq_no_DE_vector <- no_DD_MDSeq_no_DE_vector <- c()
for (one_dataset in datasets) {
  dataset_diffdist_out_dir <- sprintf("%s/%s", diffdist_out_dir, one_dataset)
  diffdist_output_basename <- sprintf("%s_%s_vs_%s_%s_DiffDist_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diffdist_results_file <- sprintf("%s/%s.csv", dataset_diffdist_out_dir, diffdist_output_basename)
  
  if (file.exists(diffdist_results_file)) {
    diffdist_results_df <- read.csv(file=diffdist_results_file, header=TRUE, row.names=1)
    
    ## DD+
    DDp <- sum(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] < pval_threshold & diffdist_results_df[[sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]] > 0, na.rm=TRUE)
    DDp_genes <- rownames(diffdist_results_df[which(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] < pval_threshold & diffdist_results_df[[sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]] > 0),])
    DDp_MDSeq_no_DE_genes <- DDp_genes[DDp_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DDp_MDSeq_no_DE <- length(DDp_MDSeq_no_DE_genes)
    all_datasets_diffdist_DD_MDSeq_non_DE_genes[[one_dataset]][["DDp_MDSeq_no_DE"]] <- DDp_MDSeq_no_DE_genes
    
    ## DD-
    DDm <- sum(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] < pval_threshold & diffdist_results_df[[sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]] < 0, na.rm=TRUE)
    DDm_genes <- rownames(diffdist_results_df[which(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] < pval_threshold & diffdist_results_df[[sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]] < 0),])
    DDm_MDSeq_no_DE_genes <- DDm_genes[DDm_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    DDm_MDSeq_no_DE <- length(DDm_MDSeq_no_DE_genes)
    all_datasets_diffdist_DD_MDSeq_non_DE_genes[[one_dataset]][["DDm_MDSeq_no_DE"]] <- DDm_MDSeq_no_DE
    
    ## non-DD
    no_DD <- sum(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] > pval_threshold, na.rm=TRUE)
    no_DD_genes <- rownames(diffdist_results_df[which(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] > pval_threshold),])
    no_DD_MDSeq_no_DE_genes <- no_DD_genes[no_DD_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    no_DD_MDSeq_no_DE <- length(no_DD_MDSeq_no_DE_genes)
    all_datasets_diffdist_DD_MDSeq_non_DE_genes[[one_dataset]][["no_DD_MDSeq_no_DE"]] <- no_DD_MDSeq_no_DE
    
    
  } else {
    DDp_MDSeq_no_DE <- DDm_MDSeq_no_DE <- no_DD_MDSeq_no_DE <- NA
  }
  # concatenate vectors
  DDp_MDSeq_no_DE_vector <- c(DDp_MDSeq_no_DE_vector, DDp_MDSeq_no_DE)
  DDm_MDSeq_no_DE_vector <- c(DDm_MDSeq_no_DE_vector, DDm_MDSeq_no_DE)
  no_DD_MDSeq_no_DE_vector <- c(no_DD_MDSeq_no_DE_vector, no_DD_MDSeq_no_DE)
}
all_datasets_diffdist_df <- rbind(DDp_MDSeq_no_DE_vector, DDm_MDSeq_no_DE_vector, no_DD_MDSeq_no_DE_vector)
colnames(all_datasets_diffdist_df) <- datasets
rownames(all_datasets_diffdist_df) <- c("DD+ MDSeq non-DE", "DD- MDSeq non-DE", "non-DD MDSeq non-DE")
write.csv(all_datasets_diffdist_df, file=sprintf("%s/%s_DiffDist_DD_MDSeq_non_DE.csv", output_dir, fig_output_basename), quote=FALSE, row.names=TRUE) 




# DE gene and DD gene counts
pdf(file=sprintf("%s/%s.pdf", output_dir, fig_output_basename), width=9)
## DE+, DE-, non-DE, total
mdseq_categories <- c("DE+", "DE-", "non-DE")
MDSeq_DE_df2ggplot <- data.frame()
for (one_dataset in colnames(all_datasets_mdseq_df)) {
  one_df <- data.frame(dataset=rep(one_dataset, length(mdseq_categories)), category=mdseq_categories, count=as.numeric(all_datasets_mdseq_df[mdseq_categories, one_dataset]), method=rep("MDSeq", length(mdseq_categories)))
  MDSeq_DE_df2ggplot <- rbind(MDSeq_DE_df2ggplot, one_df)
}
### reorder levels
MDSeq_DE_df2ggplot$category <- factor(MDSeq_DE_df2ggplot$category, levels=rev(mdseq_categories))
MDSeq_DE_plot <- ggplot(data = MDSeq_DE_df2ggplot, aes(x=dataset, y=category, fill=count)) + 
  geom_tile(alpha=0.8) +
  geom_text(aes(dataset, category, label=count), color="black", size=3) +
  scale_fill_distiller(palette = "YlOrRd", direction=1) +
  labs(fill = "Counts") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(panel.border = element_blank()) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.title=element_text(size=10)) +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=8)) +
  theme(axis.ticks = element_blank())

## DD+ MDSeq non-DE, DD- MDSeq non-DE, non-DD MDSeq non-DE
categories_counts <- function(count_df, method, method_cat, common_cat) {
  df2return <- data.frame()
  for (one_dataset in colnames(count_df)) {
    one_df <- data.frame(dataset=rep(one_dataset, length(method_cat)), category=method_cat, count=as.numeric(count_df[method_cat,one_dataset]), method=rep(method, length(method_cat)))
    ### rename categories
    one_df$category <- factor(common_cat)
    df2return <- rbind(df2return, one_df)
  }
  return(df2return)
}

Levene_categories <- c("DV+ MDSeq non-DE", "DV- MDSeq non-DE", "non-DV MDSeq non-DE")
diphiseq_categories <- gamlss_categories <- diffdist_categories <- c("DD+ MDSeq non-DE", "DD- MDSeq non-DE", "non-DD MDSeq non-DE")
mdseq_categories <- c("DV+ non-DE", "DV- non-DE", "non-DV non-DE")
common_categories <- c("DD+/DV+", "DD-/DV-", "non-DD/non-DV")
df2ggplot_barplot <- data.frame()
df2ggplot_barplot <- rbind(df2ggplot_barplot, categories_counts(all_datasets_Levene_df, "Levene", Levene_categories, common_categories))
df2ggplot_barplot <- rbind(df2ggplot_barplot, categories_counts(all_datasets_mdseq_df, "MDSeq", mdseq_categories, common_categories))
df2ggplot_barplot <- rbind(df2ggplot_barplot, categories_counts(all_datasets_diphiseq_df, "DiPhiSeq", diphiseq_categories, common_categories))
df2ggplot_barplot <- rbind(df2ggplot_barplot, categories_counts(all_datasets_gamlss_df, "GAMLSS", gamlss_categories, common_categories))
df2ggplot_barplot <- rbind(df2ggplot_barplot, categories_counts(all_datasets_diffdist_df, "DiffDist", diffdist_categories, common_categories))
### reorder levels
df2ggplot_barplot$category <- factor(df2ggplot_barplot$category, levels=common_categories)
df2ggplot_barplot$method <- factor(df2ggplot_barplot$method, levels=method_order[method_order %in% df2ggplot_barplot$method])
count_barplot <- ggplot(df2ggplot_barplot[which(df2ggplot_barplot$category %in% c("DD+/DV+", "DD-/DV-")),], aes(x=dataset, y=count, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~category) +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Dataset", y="#genes", fill="Method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12)) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10)) +
  theme(axis.text.x=element_text(angle=30, hjust=1))

### multiplot
legend_max <- max(MDSeq_DE_df2ggplot$count, na.rm=TRUE)
MDSeq_DE_df2ggplot$category <- factor(MDSeq_DE_df2ggplot$category, levels=c("non-DE", "DE-", "DE+"))
DE_DD_grid_plot <- plot_grid(MDSeq_DE_plot, count_barplot, ncol=1, align="v", labels=c("A", "B"))
print(DE_DD_grid_plot)
dev.off()


############
# Figure 5 #
############
fig_output_basename <- sprintf("fig5_S2file-TCGA_%s_vs_%s_FC_%d_upset", pheno_1, pheno_2, DD_FC_threshold)
all_datasets_upset_df <- data.frame()
figure_plot_list <- list()
supp_figure_plot_list <- list()
for (one_dataset in datasets) {
  # load data
  ## read Levene's test output table
  dataset_Levene_out_dir <- sprintf("%s/10-Levene/%s", output_dir, one_dataset)
  Levene_output_basename <- sprintf("%s_%s_vs_%s_%s_Levene_test_results", one_dataset, pheno_1, pheno_2, normalization_method)
  Levene_results_file <- sprintf("%s/%s.csv", dataset_Levene_out_dir, Levene_output_basename)
  Levene_results_df <- read.csv(Levene_results_file, header=TRUE, row.names=1)
  
  ## read MDSeq output table
  dataset_mdseq_out_dir <- sprintf("%s/20-MDSeq/%s", output_dir, one_dataset)
  mdseq_output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq_FC_%s", one_dataset, pheno_1, pheno_2, normalization_method, sub("[.]", "_", log(DD_FC_threshold, 2)))
  mdseq_results_file <- sprintf("%s/%s.csv", dataset_mdseq_out_dir, mdseq_output_basename)
  mdseq_results_df <- read.csv(mdseq_results_file, header=TRUE, row.names=1)
  
  ## read DiPhiSeq output table
  dataset_diphiseq_out_dir <- sprintf("%s/30-DiPhiSeq/%s", output_dir, one_dataset)
  diphiseq_output_basename <- sprintf("%s_%s_vs_%s_%s_DiPhiSeq_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diphiseq_results_file <- sprintf("%s/%s.csv", dataset_diphiseq_out_dir, diphiseq_output_basename)
  diphiseq_results_df <- read.csv(diphiseq_results_file, header=TRUE, row.names=1)
  
  ## read GAMLSS output table
  dataset_gamlss_out_dir <- sprintf("%s/40-GAMLSS/%s", output_dir, one_dataset)
  gamlss_output_basename <- sprintf("%s_%s_vs_%s_%s_GAMLSS_results", one_dataset, pheno_1, pheno_2, normalization_method)
  gamlss_results_file <- sprintf("%s/%s.csv", dataset_gamlss_out_dir, gamlss_output_basename)
  gamlss_results_df <- read.csv(gamlss_results_file, header=TRUE, row.names=1)
  
  ## read DiffDist output table
  dataset_diffdist_out_dir <- sprintf("%s/50-DiffDist/%s", output_dir, one_dataset)
  diffdist_output_basename <- sprintf("%s_%s_vs_%s_%s_DiffDist_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diffdist_results_file <- sprintf("%s/%s.csv", dataset_diffdist_out_dir, diffdist_output_basename)
  diffdist_results_df <- read.csv(diffdist_results_file, header=TRUE, row.names=1)
  
  # MDSeq non-DE genes
  mdseq_non_DE <- rownames(mdseq_results_df[which(mdseq_results_df$FDR.mean > 0.05),])
  ## Levene's test
  mdseq_non_DE_Levene_df <- Levene_results_df[mdseq_non_DE,]
  mdseq_non_DE_Levene_DDp <- rownames(mdseq_non_DE_Levene_df[which(mdseq_non_DE_Levene_df$p_value.BH < 0.05 & mdseq_non_DE_Levene_df[, sprintf("var_%s", pheno_1)] > mdseq_non_DE_Levene_df[, sprintf("var_%s", pheno_2)]),])
  df2ggplot <- data.frame(gene=mdseq_non_DE_Levene_DDp, method=rep("Levene", length(mdseq_non_DE_Levene_DDp)))
  ## MDSeq
  mdseq_non_DE_df <- mdseq_results_df[mdseq_non_DE,]
  mdseq_non_DE_DDp <- rownames(mdseq_non_DE_df[which(mdseq_non_DE_df$FDR.dispersion < 0.05 & mdseq_non_DE_df[, sprintf("%svs%s.dispersion.log2FC.%s", pheno_1, pheno_2, log(DD_FC_threshold, 2))] > 0),])
  df2ggplot <- rbind(df2ggplot, data.frame(gene=mdseq_non_DE_DDp, method=rep("MDSeq", length(mdseq_non_DE_DDp))))
  ## DiPhiSeq
  mdseq_non_DE_diphiseq_df <- diphiseq_results_df[mdseq_non_DE,]
  mdseq_non_DE_diphiseq_DDp <- rownames(mdseq_non_DE_diphiseq_df[which(mdseq_non_DE_diphiseq_df$fdr.phi < 0.05 & mdseq_non_DE_diphiseq_df[, sprintf("phi.%s", pheno_1)] > mdseq_non_DE_diphiseq_df[, sprintf("phi.%s", pheno_2)]),])
  df2ggplot <- rbind(df2ggplot, data.frame(gene=mdseq_non_DE_diphiseq_DDp, method=rep("DiPhiSeq", length(mdseq_non_DE_diphiseq_DDp))))
  ## GAMLSS
  mdseq_non_DE_gamlss_df <- gamlss_results_df[mdseq_non_DE,]
  mdseq_non_DE_gamlss_DDp <- rownames(mdseq_non_DE_gamlss_df[which(mdseq_non_DE_gamlss_df$padj.cv < 0.05 & mdseq_non_DE_gamlss_df[, sprintf("CV.%s", pheno_1)] > mdseq_non_DE_gamlss_df[, sprintf("CV.%s", pheno_2)]),])
  df2ggplot <- rbind(df2ggplot, data.frame(gene=mdseq_non_DE_gamlss_DDp, method=rep("GAMLSS", length(mdseq_non_DE_gamlss_DDp))))
  ## DiffDist
  mdseq_non_DE_diffdist_df <- diffdist_results_df[mdseq_non_DE,]
  mdseq_non_DE_diffdist_DDp <- rownames(mdseq_non_DE_diffdist_df[which(mdseq_non_DE_diffdist_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] < 0.05 & mdseq_non_DE_diffdist_df[[sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]] > 0),])
  df2ggplot <- rbind(df2ggplot, data.frame(gene=mdseq_non_DE_diffdist_DDp, method=rep("DiffDist", length(mdseq_non_DE_diffdist_DDp))))
  
  # prepare data frame for upset plot
  df2ggplot$gene <- as.character(df2ggplot$gene)
  df2ggplot$method <- as.character(df2ggplot$method)
  df2ggplot_2 <- df2ggplot %>%
    group_by(gene) %>%
    summarise(method = list(method))
  df2ggplot_3 <- df2ggplot_2 %>%
    group_by(gene) %>%
    summarise(method = unlist(lapply(method, function(x) {paste(x, collapse="&")})))
  all_datasets_upset_df <- rbind(all_datasets_upset_df, df2ggplot_3)
  expression_upset_input <- as.numeric(table(df2ggplot_3$method))
  names(expression_upset_input) <- names(table(df2ggplot_3$method))
  p_supp_figure <- upset(fromExpression(expression_upset_input), order.by = "freq", point.size=4, line.size=1.5)
  supp_figure_plot_list[[one_dataset]] <- p_supp_figure
  p_figure <- upset(fromExpression(expression_upset_input), order.by = "freq", point.size=2, line.size=1, mb.ratio=c(0.6, 0.4))
  figure_plot_list[[one_dataset]] <- p_figure
}
write.csv(all_datasets_upset_df, file=sprintf("%s/%s.csv", output_dir, fig_output_basename), quote=FALSE, row.names=FALSE)

fig_output_basename <- sprintf("fig5-TCGA_%s_vs_%s_FC_%d_upset", pheno_1, pheno_2, DD_FC_threshold)
pdf(sprintf("%s/%s.pdf", output_dir, fig_output_basename), width=9)
# convert UpSetR plot to ggplot plot: https://github.com/hms-dbmi/UpSetR/issues/105
figure_datasets <- c("TCGA-KIRC", "TCGA-KIRP")
dataset_1 <- figure_datasets[1]
upset_1 <- figure_plot_list[[dataset_1]]
upset_1_ggplot <- cowplot::plot_grid(NULL, upset_1$Main_bar, upset_1$Sizes, upset_1$Matrix,
                                     nrow=2, align='hv', rel_heights = c(3,1),
                                     rel_widths = c(1,3))
dataset_2 <- figure_datasets[2]
upset_2 <- figure_plot_list[[dataset_2]]
upset_2_ggplot <- cowplot::plot_grid(NULL, upset_2$Main_bar, upset_2$Sizes, upset_2$Matrix,
                                     nrow=2, align='hv', rel_heights = c(3,1),
                                     rel_widths = c(1,3))

figure_plot <- ggdraw() +
  draw_plot(upset_1_ggplot, 0, 0.5, 1, 0.5) +
  draw_plot(upset_2_ggplot, 0, 0, 1, 0.5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size=15)
print(figure_plot)
grid.text(sprintf("%s", dataset_1), x = 0.65, y=0.95, gp=gpar(fontsize=12))
grid.text(sprintf("%s", dataset_2), x = 0.65, y=0.45, gp=gpar(fontsize=12))
dev.off()

fig_output_basename <- sprintf("figS2-TCGA_%s_vs_%s_FC_%d_upset", pheno_1, pheno_2, DD_FC_threshold)
pdf(sprintf("%s/%s.pdf", output_dir, fig_output_basename), width=9)
for (one_dataset in c(datasets[! datasets %in% figure_datasets])) {
  print(supp_figure_plot_list[[one_dataset]])
  grid.text(sprintf("%s", one_dataset), x = 0.65, y=0.95, gp=gpar(fontsize=14))
}
dev.off()

