library(ggplot2)
library(RColorBrewer)
library(cowplot)


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


#############
# Functions #
#############
build_DE_DD_count_heatmap_data_frame <- function(diphiseq_df, mdseq_df, sub_categories, main_categories) {
  df2ggplot <- data.frame()
  for (one_dataset in colnames(diphiseq_df)) {
    one_df <- data.frame(dataset=rep(one_dataset, length(sub_categories)), main_category=main_categories, sub_category=sub_categories, count=as.numeric(diphiseq_df[sub_categories,one_dataset]), method=rep("DiPhiSeq", length(sub_categories)))
    df2ggplot <- rbind(df2ggplot, one_df)
  }
  for (one_dataset in colnames(mdseq_df)) {
    one_df <- data.frame(dataset=rep(one_dataset, length(sub_categories)), main_category=main_categories, sub_category=sub_categories, count=as.numeric(mdseq_df[sub_categories,one_dataset]), method=rep("MDSeq", length(sub_categories)))
    df2ggplot <- rbind(df2ggplot, one_df)
  }
  # reorder levels
  df2ggplot$sub_category <- factor(df2ggplot$sub_category, levels=sub_categories)
  return(df2ggplot)
}

DE_DD_count_heatmap <- function(df, max_value, plot_legend) {
  p <- ggplot(data = df, aes(x=dataset, y=sub_category, fill=count)) + 
    facet_grid(~method) +
    geom_tile(alpha=0.8) +
    geom_text(aes(dataset, sub_category, label = count), color = "black", size = 2.5) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=8)) +
    theme(axis.ticks = element_blank())
  if (! is.na(max_value)) {
    p <- p + scale_fill_distiller(palette = "YlOrRd", direction=1, limits=c(0, max_value))
  } else {
    p <- p + scale_fill_distiller(palette = "YlOrRd", direction=1)
  }
  if (plot_legend) {
    p <- p + theme(legend.position = "bottom") +
      labs(fill = "Counts") +
      theme(legend.title=element_text(size=10)) +
      theme(legend.text=element_text(size=8))
  } else {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

build_common_specific_barplot_data_frame <- function(MDSeq_gene_list, DiPhiSeq_gene_list, gene_category) {
  dataset_vector <- count_vector <- category_vector <- percentage_vector <- c()
  datasets <- names(MDSeq_gene_list)
  for (one_dataset in datasets) {
    mdseq_genes <- MDSeq_gene_list[[one_dataset]][[gene_category]]
    diphiseq_genes <- DiPhiSeq_gene_list[[one_dataset]][[gene_category]]
    # common genes
    common_genes <- mdseq_genes[mdseq_genes %in% diphiseq_genes]
    common_genes_nb <- length(common_genes)
    # DiPhiSeq specific genes
    diphiseq_specific_genes <- diphiseq_genes[! diphiseq_genes %in% common_genes]
    diphiseq_specific_genes_nb <- length(diphiseq_specific_genes)
    # MDSeq specific genes
    mdseq_specific_genes <- mdseq_genes[! mdseq_genes %in% common_genes]
    mdseq_specific_genes_nb <- length(mdseq_specific_genes)
    # total genes
    total_genes <- common_genes_nb + diphiseq_specific_genes_nb + mdseq_specific_genes_nb
    
    dataset_vector <- c(dataset_vector, rep(one_dataset, 3))
    count_vector <- c(count_vector, c(common_genes_nb, diphiseq_specific_genes_nb, mdseq_specific_genes_nb))
    category_vector <- c(category_vector, c("DiPhiSeq\nand MDSeq", "DiPhiSeq only", "MDSeq only"))
    percentage_vector <- c(percentage_vector, c(common_genes_nb/total_genes*100, diphiseq_specific_genes_nb/total_genes*100, mdseq_specific_genes_nb/total_genes*100))
  }
  df2return <- data.frame(dataset=dataset_vector, count=count_vector, pct=percentage_vector, category=category_vector)
  df2return$category <- factor(df2return$category, levels=c("MDSeq only", "DiPhiSeq\nand MDSeq", "DiPhiSeq only"))
  df2return$dataset <- factor(df2return$dataset, levels=rev(sort(levels(df2return$dataset))), ordered=TRUE)
  return(df2return)
}

common_specific_barplot <- function(df) {
  p <- ggplot(df, aes(x=dataset, y=count, fill=category)) +
    geom_bar(stat="identity") +
    coord_flip() +
    geom_text(aes(label=sprintf("%.1f%%",pct)), size=4, position=position_stack(vjust=0.5)) +
    scale_fill_manual(values=c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[4], brewer.pal(9, "Set1")[1])) +
    labs(y="# genes") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_blank(), legend.title=element_blank()) +
    theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=10, hjust=1.5), legend.text=element_text(size=11)) +
    theme(axis.ticks.y=element_blank()) +
    guides(fill=guide_legend(reverse=TRUE))
  return(p)
}


############
# Analysis #
############
all_genes_output_basename <- sprintf("fig4_S4_S5-TCGA_%s_vs_%s_FC_%d_DE_DD_counts", pheno_1, pheno_2, DD_FC_threshold)
lowly_DE_genes_output_basename <- sprintf("figS6_S7-TCGA_%s_vs_%s_lowly_DE_DD_counts", pheno_1, pheno_2)

# MDSeq results
mdseq_out_dir <- sprintf("%s/20-MDSeq", output_dir)
## all genes
DE_vector <- DEp_vector <- DEm_vector <- DEp_DDp_vector <- DEp_no_DD_vector <- DEp_DDm_vector <- DEm_DDp_vector <- DEm_no_DD_vector <- DEm_DDm_vector <- c()
DD_vector <- DDp_vector <- DDm_vector <- DDp_no_DE_vector <- DDm_no_DE_vector <- c()
no_DD_vector <- no_DE_vector <- no_DE_no_DD_vector <- c()
all_genes_total_vector <- c()
## lowly DE genes
## count genes according differential expression and differential dispersion categories and get lowly DE genes
highly_DE_vector <- highly_DEp_vector <- highly_DEm_vector <- lowly_DE_vector <- c()
DDp_lowly_DE_vector <- DDm_lowly_DE_vector <- lowly_DE_no_DD_vector <- c()
lowly_DE_genes_total_vector <- c()
all_datasets_MDSeq_genes <- list()
for (one_dataset in datasets) {
  dataset_mdseq_out_dir <- sprintf("%s/%s", mdseq_out_dir, one_dataset)
  ## all gene analysis
  ### read MDSeq output table
  dataset_mdseq_all_genes_out_dir <- sprintf("%s/00-all_genes", dataset_mdseq_out_dir)
  mdseq_output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq_FC_%s", one_dataset, pheno_1, pheno_2, normalization_method, sub("[.]", "_", DD_FC_threshold))
  mdseq_all_genes_file <- sprintf("%s/%s.csv", dataset_mdseq_all_genes_out_dir, mdseq_output_basename)
  if (file.exists(mdseq_all_genes_file)) {
    mdseq_all_genes_df <- read.csv(file=mdseq_all_genes_file, header=TRUE, row.names=1)
    ### count genes
    #### total
    nb_total <- dim(mdseq_all_genes_df)[1]
    #### DE
    nb_DE <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold, na.rm=TRUE)
    ##### DE+
    DEp <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ###### DE+ DD+
    DEp_DDp <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0 & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ###### DE+ non-DD
    DEp_no_DD <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df$FDR.dispersion > pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ###### DE+ DD-
    DEp_DDm <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0 & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    ##### DE-
    DEm <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    ###### DE- DD+
    DEm_DDp <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0 & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ###### DE- non-DD
    DEm_no_DD <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df$FDR.dispersion > pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    ###### DE- DD-
    DEm_DDm <- sum(mdseq_all_genes_df$FDR.mean < pval_threshold & mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.mean.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0 & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    #### DD
    nb_DD <- sum(mdseq_all_genes_df$FDR.dispersion < pval_threshold, na.rm=TRUE)
    ##### DD+
    DDp <- sum(mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ###### DD+ non-DE
    DDp_no_DE <- sum(mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df$FDR.mean > pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ##### DD-
    DDm <- sum(mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    ###### DD- non-DE
    DDm_no_DE <- sum(mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df$FDR.mean > pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    #### non-DD
    no_DD <- sum(mdseq_all_genes_df$FDR.dispersion > pval_threshold, na.rm=TRUE)
    #### non-DE
    no_DE <- sum(mdseq_all_genes_df$FDR.mean > pval_threshold, na.rm=TRUE)
    ##### non-DE and non-DD
    no_DE_no_DD <- sum(mdseq_all_genes_df$FDR.mean > pval_threshold & mdseq_all_genes_df$FDR.dispersion > pval_threshold, na.rm=TRUE)
    
    ### get genes
    #### non-DE genes
    non_DE_genes <- rownames(mdseq_all_genes_df[which(mdseq_all_genes_df$FDR.mean > pval_threshold),])
    all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]] <- non_DE_genes
    #### DD+ genes
    DDp_genes <- rownames(mdseq_all_genes_df[which(mdseq_all_genes_df$FDR.dispersion < pval_threshold & mdseq_all_genes_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0),])
    ##### DD+ among non-DE genes
    DDp_non_DE_genes <- DDp_genes[DDp_genes %in% non_DE_genes]
    all_datasets_MDSeq_genes[[one_dataset]][["DDp_non_DE"]] <- DDp_non_DE_genes
  } else {
    nb_DE <- DEp <- DEp_DDp <- DEp_no_DD <- DEp_DDm <- DEm <- DEm_DDp <- DEm_no_DD <- DEm_DDm <-NA
    nb_DD <- DDp <- DDp_no_DE <- DDm <- DDm_no_DE <- NA
    no_DD <- no_DE <- no_DE_no_DD <- nb_total <- NA
  }
  DE_vector <- c(DE_vector, nb_DE)
  DEp_vector <- c(DEp_vector, DEp)
  DEp_DDp_vector <- c(DEp_DDp_vector, DEp_DDp)
  DEp_no_DD_vector <- c(DEp_no_DD_vector, DEp_no_DD)
  DEp_DDm_vector <- c(DEp_DDm_vector, DEp_DDm)
  DEm_vector <- c(DEm_vector, DEm)
  DEm_DDp_vector <- c(DEm_DDp_vector, DEm_DDp)
  DEm_no_DD_vector <- c(DEm_no_DD_vector, DEm_no_DD)
  DEm_DDm_vector <- c(DEm_DDm_vector, DEm_DDm)
  DD_vector <- c(DD_vector, nb_DD)
  DDp_vector <- c(DDp_vector, DDp)
  DDp_no_DE_vector <- c(DDp_no_DE_vector, DDp_no_DE)
  DDm_vector <- c(DDm_vector, DDm)
  DDm_no_DE_vector <- c(DDm_no_DE_vector, DDm_no_DE)
  no_DD_vector <- c(no_DD_vector, no_DD)
  no_DE_vector <- c(no_DE_vector, no_DE)
  no_DE_no_DD_vector <- c(no_DE_no_DD_vector, no_DE_no_DD)
  all_genes_total_vector <- c(all_genes_total_vector, nb_total)
  
  ## lowly DE gene analysis
  ### fold-change threshold to identify highly genes (not in log2 scale)
  DE_FC_threshold <- ifelse(one_dataset=="TCGA-BRCA", 1.25, 1.3)
  ### highly DE and lowly DE genes
  #### read MDSeq output table
  dataset_mdseq_lowly_DE_out_dir <- sprintf("%s/10-lowly_DE_genes", dataset_mdseq_out_dir)
  mdseq_output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq_FC_%s", one_dataset, pheno_1, pheno_2, normalization_method, sub("[.]", "_", DE_FC_threshold))
  mdseq_highly_DE_file <- sprintf("%s/%s.csv", dataset_mdseq_lowly_DE_out_dir, mdseq_output_basename)
  if (file.exists(mdseq_highly_DE_file)) {
    mdseq_highly_de_df <- read.csv(file=mdseq_highly_DE_file, header=TRUE, row.names=1)
    #### count genes
    ##### total
    nb_total <- dim(mdseq_highly_de_df)[1]
    ##### highly DE
    nb_highly_DE <- sum(mdseq_highly_de_df$FDR.mean < pval_threshold, na.rm=TRUE)
    ###### highly DE+
    highly_DEp <- sum(mdseq_highly_de_df$FDR.mean < pval_threshold & mdseq_highly_de_df[[sprintf("%svs%s.mean.log2FC.%s", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ###### highly DE-
    highly_DEm <- sum(mdseq_highly_de_df$FDR.mean < pval_threshold & mdseq_highly_de_df[[sprintf("%svs%s.mean.log2FC.%s", make.names(pheno_1), make.names(pheno_2), log(DE_FC_threshold, 2))]] < 0, na.rm=TRUE)
    ##### lowly DE
    lowly_DE <- sum(mdseq_highly_de_df$FDR.mean > pval_threshold, na.rm=TRUE)
    
    #### get genes
    ##### lowly DE genes
    lowly_DE_genes <- rownames(mdseq_highly_de_df[which(mdseq_highly_de_df$FDR.mean > pval_threshold),])
    all_datasets_MDSeq_genes[[one_dataset]][["lowly_DE"]] <- lowly_DE_genes
  } else {
    nb_highly_DE <- highly_DEp <- highly_DEm <- lowly_DE <- nb_total <-NA
  }
  highly_DE_vector <- c(highly_DE_vector, nb_highly_DE)
  highly_DEp_vector <- c(highly_DEp_vector, highly_DEp)
  highly_DEm_vector <- c(highly_DEm_vector, highly_DEm)
  lowly_DE_vector <- c(lowly_DE_vector, lowly_DE)
  lowly_DE_genes_total_vector <- c(lowly_DE_genes_total_vector, nb_total) 
  
  ### DD+ genes among lowly DE genes
  #### read MDSeq output table
  mdseq_output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq_lowly_DE_FC_%s", one_dataset, pheno_1, pheno_2, normalization_method, sub("[.]", "_", DD_FC_threshold))
  mdseq_lowly_DE_file <- sprintf("%s/%s.csv", dataset_mdseq_lowly_DE_out_dir, mdseq_output_basename)
  if (file.exists(mdseq_lowly_DE_file)) {
    mdseq_lowly_de_df <- read.csv(file=mdseq_lowly_DE_file, header=TRUE, row.names=1)
    #### count genes
    ##### DD+ lowly DE
    DDp_lowly_DE <- sum(mdseq_lowly_de_df$FDR.dispersion < pval_threshold & mdseq_lowly_de_df[[sprintf("%svs%s.dispersion.log2FC.%s", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0, na.rm=TRUE)
    ##### DD- lowly DE
    DDm_lowly_DE <- sum(mdseq_lowly_de_df$FDR.dispersion < pval_threshold & mdseq_lowly_de_df[[sprintf("%svs%s.dispersion.log2FC.%s", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] < 0, na.rm=TRUE)
    ##### no DD lowly DE
    lowly_DE_no_DD <- sum(mdseq_lowly_de_df$FDR.dispersion > pval_threshold, na.rm=TRUE)
    
    #### get genes
    ##### DD+ genes among lowly DE genes
    DDp_lowly_DE_genes <- rownames(mdseq_lowly_de_df[which(mdseq_lowly_de_df$FDR.dispersion < pval_threshold & mdseq_lowly_de_df[[sprintf("%svs%s.dispersion.log2FC.%d", make.names(pheno_1), make.names(pheno_2), log(DD_FC_threshold, 2))]] > 0),])
    all_datasets_MDSeq_genes[[one_dataset]][["DDp_lowly_DE"]] <- DDp_lowly_DE_genes
  } else {
    nb_DD <- DDp_lowly_DE <- DDm_lowly_DE <- lowly_DE_no_DD <- NA
  }
  DDp_lowly_DE_vector <- c(DDp_lowly_DE_vector, DDp_lowly_DE)
  DDm_lowly_DE_vector <- c(DDm_lowly_DE_vector, DDm_lowly_DE)
  lowly_DE_no_DD_vector <- c(lowly_DE_no_DD_vector, lowly_DE_no_DD)
}
all_datasets_mdseq_all_genes_df <- rbind(DE_vector, DEp_vector, DEp_DDp_vector, DEp_no_DD_vector, DEp_DDm_vector, DEm_vector, DEm_DDp_vector, DEm_no_DD_vector, DEm_DDm_vector, DD_vector, DDp_vector, DDp_no_DE_vector, DDm_vector, DDm_no_DE_vector, no_DD_vector, no_DE_vector, no_DE_no_DD_vector, all_genes_total_vector)
colnames(all_datasets_mdseq_all_genes_df) <- datasets
rownames(all_datasets_mdseq_all_genes_df) <- c("DE", "DE+", "DD+ DE+", "non-DD DE+", "DD- DE+", "DE-", "DD+ DE-", "non-DD DE-", "DD- DE-", "DD", "DD+", "DD+ non-DE", "DD-", "DD- non-DE", "non-DD", "non-DE", "non-DD non-DE", "total")
write.csv(all_datasets_mdseq_all_genes_df, file=sprintf("%s/%s_MDSeq.csv", mdseq_out_dir, all_genes_output_basename), quote=FALSE, row.names=TRUE) 
all_datasets_mdseq_lowly_DE_df <- rbind(highly_DE_vector, highly_DEp_vector, highly_DEm_vector, DDp_lowly_DE_vector, DDm_lowly_DE_vector, lowly_DE_vector, lowly_DE_no_DD_vector, lowly_DE_genes_total_vector)
colnames(all_datasets_mdseq_lowly_DE_df) <- datasets
rownames(all_datasets_mdseq_lowly_DE_df) <- c("highly DE", "highly DE+", "highly DE-", "DD+ lowly DE", "DD- lowly DE", "lowly DE", "no DD lowly DE", "total")
write.csv(all_datasets_mdseq_lowly_DE_df, file=sprintf("%s/%s_MDSeq.csv", mdseq_out_dir, lowly_DE_genes_output_basename), quote=FALSE, row.names=TRUE) 

# DiPhiSeq results
DD_vector <- DDp_vector <- DDm_vector <- DDp_DEp_vector <- DDp_no_DE_vector <- DDp_DEm_vector <- DDm_DEp_vector <- DDm_no_DE_vector <- DDm_DEm_vector <- c()
DE_vector <- DEp_vector <- DEm_vector <- no_DD_DEp_vector <- no_DD_DEm_vector <- c()
no_DD_vector <- no_DE_vector <- no_DD_no_DE_vector <- c()
total_vector <- c()
## count genes according differential dispersion among lowly DE genes identified by MDSeq
DDp_MDSeq_lowly_DE_vector <- DDm_MDSeq_lowly_DE_vector <- no_DD_MDSeq_lowly_DE_vector <- c()
all_datasets_DiPhiSeq_DDp_genes <- list()
for (one_dataset in datasets) {
  ## read DiPhiSeq output table
  diphiseq_out_dir <- sprintf("%s/10-DiPhiSeq/%s", output_dir, one_dataset)
  diphiseq_output_basename <- sprintf("%s_%s_vs_%s_%s_DiPhiSeq_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diphiseq_file <- sprintf("%s/%s.csv", diphiseq_out_dir, diphiseq_output_basename)
  if (file.exists(diphiseq_file)) {
    diphiseq_df <- read.csv(file=diphiseq_file, header=TRUE, row.names=1)
    ## count genes
    ### total
    total <- dim(diphiseq_df)[1]
    ### DD
    DD <- sum(diphiseq_df$fdr.phi < pval_threshold, na.rm=TRUE)
    #### DD+
    DDp <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] > diphiseq_df[, sprintf("phi.%s", pheno_2)], na.rm=TRUE)
    ##### DD+ DE+
    DDp_DEp <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] > diphiseq_df[, sprintf("phi.%s", pheno_2)] & diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] > diphiseq_df[, sprintf("beta.%s", pheno_2)], na.rm=TRUE)
    ##### DD+ non-DE
    DDp_no_DE <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] > diphiseq_df[, sprintf("phi.%s", pheno_2)] & diphiseq_df$fdr.beta > pval_threshold, na.rm=TRUE)
    ##### DD+ DE-
    DDp_DEm <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] > diphiseq_df[, sprintf("phi.%s", pheno_2)] & diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] < diphiseq_df[, sprintf("beta.%s", pheno_2)], na.rm=TRUE)
    #### DD-
    DDm <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] < diphiseq_df[, sprintf("phi.%s", pheno_2)], na.rm=TRUE)
    ##### DD- DE+
    DDm_DEp <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] < diphiseq_df[, sprintf("phi.%s", pheno_2)] & diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] > diphiseq_df[, sprintf("beta.%s", pheno_2)], na.rm=TRUE)
    ##### DD- non-DE
    DDm_no_DE <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] < diphiseq_df[, sprintf("phi.%s", pheno_2)] & diphiseq_df$fdr.beta > pval_threshold, na.rm=TRUE)
    ##### DD- DE-
    DDm_DEm <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] < diphiseq_df[, sprintf("phi.%s", pheno_2)] & diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] < diphiseq_df[, sprintf("beta.%s", pheno_2)], na.rm=TRUE)
    ### DE
    DE <- sum(diphiseq_df$fdr.beta < pval_threshold, na.rm=TRUE)
    #### DE+
    DEp <- sum(diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] > diphiseq_df[, sprintf("beta.%s", pheno_2)], na.rm=TRUE)
    ##### non-DD DE+
    no_DD_DEp <- sum(diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] > diphiseq_df[, sprintf("beta.%s", pheno_2)] & diphiseq_df$fdr.phi > pval_threshold, na.rm=TRUE)
    #### DE-
    DEm <- sum(diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] < diphiseq_df[, sprintf("beta.%s", pheno_2)], na.rm=TRUE)
    ##### non-DD DE-
    no_DD_DEm <- sum(diphiseq_df$fdr.beta < pval_threshold & diphiseq_df[, sprintf("beta.%s", pheno_1)] < diphiseq_df[, sprintf("beta.%s", pheno_2)] & diphiseq_df$fdr.phi > pval_threshold, na.rm=TRUE)
    ### non-DD
    no_DD <- sum(diphiseq_df$fdr.phi > pval_threshold, na.rm=TRUE)
    ### non-DE
    no_DE <- sum(diphiseq_df$fdr.beta > pval_threshold, na.rm=TRUE)
    #### non-DD non-DE
    no_DD_no_DE <- sum(diphiseq_df$fdr.phi > pval_threshold & diphiseq_df$fdr.beta > pval_threshold, na.rm=TRUE)
    
    ## get genes
    ### non-DE genes
    non_DE_genes <- rownames(diphiseq_df[which(diphiseq_df$fdr.beta > pval_threshold),])
    ### DD+ genes
    DDp_genes <- rownames(diphiseq_df[which(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] > diphiseq_df[, sprintf("phi.%s", pheno_2)]),])
    #### DD+ among non-DE genes identified by MDSeq
    DDp_MDSeq_non_DE_genes <- DDp_genes[DDp_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["non_DE"]]]
    all_datasets_DiPhiSeq_DDp_genes[[one_dataset]][["DDp_non_DE"]] <- DDp_MDSeq_non_DE_genes
    #### DD+ among lowly DE genes identified by MDSeq
    DDp_MDSeq_lowly_DE_genes <- DDp_genes[DDp_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["lowly_DE"]]]
    DDp_MDSeq_lowly_DE <- length(DDp_MDSeq_lowly_DE_genes)
    all_datasets_DiPhiSeq_DDp_genes[[one_dataset]][["DDp_lowly_DE"]] <- DDp_MDSeq_lowly_DE_genes
    ### DD- among lowly DE genes identified by MDSeq
    DDm <- sum(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] < diphiseq_df[, sprintf("phi.%s", pheno_2)], na.rm=TRUE)
    DDm_genes <- rownames(diphiseq_df[which(diphiseq_df$fdr.phi < pval_threshold & diphiseq_df[, sprintf("phi.%s", pheno_1)] < diphiseq_df[, sprintf("phi.%s", pheno_2)]),])
    DDm_MDSeq_lowly_DE_genes <- DDm_genes[DDm_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["lowly_DE"]]]
    DDm_MDSeq_lowly_DE <- length(DDm_MDSeq_lowly_DE_genes)
    ### no DD among lowly DE genes identified by MDSeq
    no_DD <- sum(diphiseq_df$fdr.phi > pval_threshold, na.rm=TRUE)
    no_DD_genes <- rownames(diphiseq_df[which(diphiseq_df$fdr.phi > pval_threshold),])
    no_DD_MDSeq_lowly_DE_genes <- no_DD_genes[no_DD_genes %in% all_datasets_MDSeq_genes[[one_dataset]][["lowly_DE"]]]
    no_DD_MDSeq_lowly_DE <- length(no_DD_MDSeq_lowly_DE_genes)
  } else {
    DD <- DDp <- DDp_DEp <- DDp_no_DE <- DDp_DEm <- DDm <- DDm_DEp <- DDm_no_DE <- DDm_DEm <-NA
    DE <- DEp <- no_DD_DEp <- DEm <- no_DD_DEm <- NA
    no_DD <- no_DE <- no_DD_no_DE <- total <- NA
    DDp_MDSeq_lowly_DE <- DDm_MDSeq_lowly_DE <- no_DD_MDSeq_lowly_DE <- NA
  }
  DD_vector <- c(DD_vector, DD)
  DDp_vector <- c(DDp_vector, DDp)
  DDp_DEp_vector <- c(DDp_DEp_vector, DDp_DEp)
  DDp_no_DE_vector <- c(DDp_no_DE_vector, DDp_no_DE)
  DDp_DEm_vector <- c(DDp_DEm_vector, DDp_DEm)
  DDm_vector <- c(DDm_vector, DDm)
  DDm_DEp_vector <- c(DDm_DEp_vector, DDm_DEp)
  DDm_no_DE_vector <- c(DDm_no_DE_vector, DDm_no_DE)
  DDm_DEm_vector <- c(DDm_DEm_vector, DDm_DEm)
  DE_vector <- c(DE_vector, DE)
  DEp_vector <- c(DEp_vector, DEp)
  no_DD_DEp_vector <- c(no_DD_DEp_vector, no_DD_DEp)
  DEm_vector <- c(DEm_vector, DEm)
  no_DD_DEm_vector <- c(no_DD_DEm_vector, no_DD_DEm)
  no_DD_vector <- c(no_DD_vector, no_DD)
  no_DE_vector <- c(no_DE_vector, no_DE)
  no_DD_no_DE_vector <- c(no_DD_no_DE_vector, no_DD_no_DE)
  total_vector <- c(total_vector, total)
  DDp_MDSeq_lowly_DE_vector <- c(DDp_MDSeq_lowly_DE_vector, DDp_MDSeq_lowly_DE)
  DDm_MDSeq_lowly_DE_vector <- c(DDm_MDSeq_lowly_DE_vector, DDm_MDSeq_lowly_DE)
  no_DD_MDSeq_lowly_DE_vector <- c(no_DD_MDSeq_lowly_DE_vector, no_DD_MDSeq_lowly_DE)
}
diphiseq_out_dir <- sprintf("%s/10-DiPhiSeq", output_dir) # TODO: remove before publication
all_datasets_diphiseq_df <- rbind(DD_vector, DDp_vector, DDp_DEp_vector, DDp_no_DE_vector, DDp_DEm_vector, DDm_vector, DDm_DEp_vector, DDm_no_DE_vector, DDm_DEm_vector, DE_vector, DEp_vector, no_DD_DEp_vector, DEm_vector, no_DD_DEm_vector, no_DD_vector, no_DE_vector, no_DD_no_DE_vector, total_vector)
colnames(all_datasets_diphiseq_df) <- datasets
rownames(all_datasets_diphiseq_df) <- c("DD", "DD+", "DD+ DE+", "DD+ non-DE", "DD+ DE-", "DD-", "DD- DE+", "DD- non-DE", "DD- DE-", "DE", "DE+", "non-DD DE+", "DE-", "non-DD DE-", "non-DD", "non-DE", "non-DD non-DE", "total")
write.csv(all_datasets_diphiseq_df, file=sprintf("%s/%s_DiPhiSeq.csv", diphiseq_out_dir, all_genes_output_basename), quote=FALSE, row.names=TRUE) 
all_datasets_diphiseq_MDSeq_lowly_DE_df <- rbind(DDp_MDSeq_lowly_DE_vector, DDm_MDSeq_lowly_DE_vector, no_DD_MDSeq_lowly_DE_vector)
colnames(all_datasets_diphiseq_MDSeq_lowly_DE_df) <- datasets
rownames(all_datasets_diphiseq_MDSeq_lowly_DE_df) <- c("DD+ lowly DE", "DD- lowly DE", "no DD lowly DE")
write.csv(all_datasets_diphiseq_MDSeq_lowly_DE_df, file=sprintf("%s/%s_DiPhiSeq.csv", diphiseq_out_dir, lowly_DE_genes_output_basename), quote=FALSE, row.names=TRUE) 

# DE DD count heatmaps
## all gene anlysis
pdf(file=sprintf("%s/%s.pdf", output_dir, all_genes_output_basename))
### DE, DE+, DE-, non-DE, DD+, DD-, non-DD
DE_DD_sub_categories <- c("DE", "DE+", "DE-", "non-DE", "DD", "DD+", "DD-", "non-DD")
DE_DD_categories <- c(rep("DE", 4), rep("DD", 4))
df2ggplot <- build_DE_DD_count_heatmap_data_frame(all_datasets_diphiseq_df, all_datasets_mdseq_all_genes_df, DE_DD_sub_categories, DE_DD_categories)
legend_max <- max(df2ggplot$count, na.rm=TRUE)
#### DE subplot
##### change category level order
df2ggplot_DE <- droplevels(df2ggplot[which(df2ggplot$main_category=="DE"),])
df2ggplot_DE$sub_category <- factor(df2ggplot_DE$sub_category, levels=c("non-DE", "DE-", "DE+", "DE"))
##### plot
DE_plot <- DE_DD_count_heatmap(df2ggplot_DE, legend_max, TRUE)
#### DD subplot
##### change category level order
df2ggplot_DD <- droplevels(df2ggplot[which(df2ggplot$main_category=="DD"),])
df2ggplot_DD$sub_category <- factor(df2ggplot_DD$sub_category, levels=c("non-DD", "DD-", "DD+", "DD"))
##### plot
DD_plot <- DE_DD_count_heatmap(df2ggplot_DD, legend_max, FALSE)
#### DE and DD plot
plot_legend <- get_legend(DE_plot)
DE_plot <- DE_plot + theme(legend.position = "none")
DE_DD_grid_plot <- plot_grid(DE_plot, DD_plot, ncol=1, align="v", labels=c("A", "B"))
multiplot_facet <- ggdraw() +
  draw_plot(DE_DD_grid_plot, 0, 0.1, 1, 0.9) +
  draw_plot(plot_legend, 0, 0, 1, 0.1)
print(multiplot_facet)

### DD, DD+, DD+ DE+, DD+ non-DE, DD+ DE-, DD-, DD- DE+, DD- non-DE, DD- DE-
DE_DD_sub_categories <- c("DD", "DD+", "DD+ DE+", "DD+ non-DE", "DD+ DE-", "DD-", "DD- DE+", "DD- non-DE", "DD- DE-")
DD_categories <- c("DD", rep("DD+", 4), rep("DD-", 4))
df2ggplot <- build_DE_DD_count_heatmap_data_frame(all_datasets_diphiseq_df, all_datasets_mdseq_all_genes_df, DE_DD_sub_categories, DD_categories)
legend_max <- max(df2ggplot$count, na.rm=TRUE)
#### DD+ subplot
##### change category level order and label
df2ggplot_DDp <- droplevels(df2ggplot[which(df2ggplot$main_category=="DD+"),])
df2ggplot_DDp$sub_category <- factor(df2ggplot_DDp$sub_category, levels=c("DD+ DE-", "DD+ non-DE", "DD+ DE+", "DD+"))
levels(df2ggplot_DDp$sub_category)[which(levels(df2ggplot_DDp$sub_category) == "DD+ DE-")] <- "DD+ &\nDE-"
levels(df2ggplot_DDp$sub_category)[which(levels(df2ggplot_DDp$sub_category) == "DD+ non-DE")] <- "DD+ &\nnon-DE"
levels(df2ggplot_DDp$sub_category)[which(levels(df2ggplot_DDp$sub_category) == "DD+ DE+")] <- "DD+ &\nDE+"
##### plot
DDp_plot <- DE_DD_count_heatmap(df2ggplot_DDp, legend_max, TRUE)
#### DD- subplot
##### change category level order and label
df2ggplot_DDm <- droplevels(df2ggplot[which(df2ggplot$main_category=="DD-"),])
df2ggplot_DDm$sub_category <- factor(df2ggplot_DDm$sub_category, levels=c("DD- DE-", "DD- non-DE", "DD- DE+", "DD-"))
levels(df2ggplot_DDm$sub_category)[which(levels(df2ggplot_DDm$sub_category) == "DD- DE-")] <- "DD- &\nDE-"
levels(df2ggplot_DDm$sub_category)[which(levels(df2ggplot_DDm$sub_category) == "DD- non-DE")] <- "DD- &\nnon-DE"
levels(df2ggplot_DDm$sub_category)[which(levels(df2ggplot_DDm$sub_category) == "DD- DE+")] <- "DD- &\nDE+"
##### plot
DDm_plot <- DE_DD_count_heatmap(df2ggplot_DDm, legend_max, FALSE)
#### DD+ and DD- plot
plot_legend <- get_legend(DDp_plot)
DDp_plot <- DDp_plot + theme(legend.position = "none")
DDp_DDm_grid_plot <- plot_grid(DDp_plot, DDm_plot, ncol=1, align="v", labels=c("A", "B"))
multiplot_facet <- ggdraw() +
  draw_plot(DDp_DDm_grid_plot, 0, 0.1, 1, 0.9) +
  draw_plot(plot_legend, 0, 0, 1, 0.1)
print(multiplot_facet)

### DE+, DE-, non-DE, DD+ non-DE, DD- non-DE, non-DD non-DE
DE_DD_sub_categories <- c("DE+", "DE-", "non-DE", "DD+ non-DE", "DD- non-DE", "non-DD non-DE")
DE_DD_categories <- c(rep("DE", 3), rep("DD", 3))
df2ggplot <- build_DE_DD_count_heatmap_data_frame(all_datasets_diphiseq_df, all_datasets_mdseq_all_genes_df, DE_DD_sub_categories, DE_DD_categories)
legend_max <- max(df2ggplot$count, na.rm=TRUE)
#### DE subplot
##### change category level order
df2ggplot_DE <- droplevels(df2ggplot[which(df2ggplot$main_category=="DE"),])
df2ggplot_DE$sub_category <- factor(df2ggplot_DE$sub_category, levels=c("non-DE", "DE-", "DE+"))
##### plot
DE_plot <- DE_DD_count_heatmap(df2ggplot_DE, legend_max, TRUE)
#### DD non DE subplot
##### change category level order and label
df2ggplot_DD <- droplevels(df2ggplot[which(df2ggplot$main_category=="DD"),])
df2ggplot_DD$sub_category <- factor(df2ggplot_DD$sub_category, levels=c("non-DD non-DE", "DD- non-DE", "DD+ non-DE"))
levels(df2ggplot_DD$sub_category)[which(levels(df2ggplot_DD$sub_category) == "non-DD non-DE")] <- "non-DD &\nnon-DE"
levels(df2ggplot_DD$sub_category)[which(levels(df2ggplot_DD$sub_category) == "DD- non-DE")] <- "DD- &\nnon-DE"
levels(df2ggplot_DD$sub_category)[which(levels(df2ggplot_DD$sub_category) == "DD+ non-DE")] <- "DD+ &\nnon-DE"
##### plot
DD_plot <- DE_DD_count_heatmap(df2ggplot_DD, legend_max, FALSE)
#### DE and DD non DE plot
plot_legend <- get_legend(DE_plot)
DE_plot <- DE_plot + theme(legend.position = "none")
DE_DD_grid_plot <- plot_grid(DE_plot, DD_plot, ncol=1, align="v", labels=c("A", "B"))
multiplot_facet <- ggdraw() +
  draw_plot(DE_DD_grid_plot, 0, 0.1, 1, 0.9) +
  draw_plot(plot_legend, 0, 0, 1, 0.1)
print(multiplot_facet)
dev.off()

## lowly DE gene analysis
lowly_DE_genes_output_basename <- sprintf("figS6-TCGA_%s_vs_%s_lowly_DE_DD_counts", pheno_1, pheno_2)
pdf(file=sprintf("%s/%s.pdf", output_dir, lowly_DE_genes_output_basename), width=9)
### highly DE and lowly DE subplot
mdseq_categories <- c("highly DE+", "highly DE-", "lowly DE")
MDSeq_DE_df2ggplot <- data.frame()
for (one_dataset in colnames(all_datasets_mdseq_lowly_DE_df)) {
  one_df <- data.frame(dataset=rep(one_dataset, length(mdseq_categories)), sub_category=mdseq_categories, count=as.numeric(all_datasets_mdseq_lowly_DE_df[mdseq_categories,one_dataset]), method=rep("MDSeq", length(mdseq_categories)))
  MDSeq_DE_df2ggplot <- rbind(MDSeq_DE_df2ggplot, one_df)
}
#### reorder levels
MDSeq_DE_df2ggplot$sub_category <- factor(MDSeq_DE_df2ggplot$sub_category, levels=c("lowly DE", "highly DE-", "highly DE+"))
legend_max <- max(MDSeq_DE_df2ggplot$count, na.rm=TRUE)
#### plot
DE_plot <- DE_DD_count_heatmap(MDSeq_DE_df2ggplot, legend_max, TRUE)
### DD lowly DE subplot
DD_sub_categories <- c("DD+ lowly DE", "DD- lowly DE", "no DD lowly DE")
DD_categories <- c("DD+", "DD-", "no DD")
df2ggplot_DD <- build_DE_DD_count_heatmap_data_frame(all_datasets_diphiseq_MDSeq_lowly_DE_df, all_datasets_mdseq_lowly_DE_df, DD_sub_categories, DD_categories)
#### change category level order and label
df2ggplot_DD$sub_category <- factor(df2ggplot_DD$sub_category, levels=c("no DD lowly DE", "DD- lowly DE", "DD+ lowly DE"))
levels(df2ggplot_DD$sub_category)[which(levels(df2ggplot_DD$sub_category) == "no DD lowly DE")] <- "no DD &\nlowly DE"
levels(df2ggplot_DD$sub_category)[which(levels(df2ggplot_DD$sub_category) == "DD- lowly DE")] <- "DD- &\nlowly DE"
levels(df2ggplot_DD$sub_category)[which(levels(df2ggplot_DD$sub_category) == "DD+ lowly DE")] <- "DD+ &\nlowly DE"
#### plot
DD_plot <- DE_DD_count_heatmap(df2ggplot_DD, legend_max, FALSE)
### highly DE, lowly DE and DD lowly DE plot
plot_legend <- get_legend(DE_plot)
DE_plot <- DE_plot + theme(legend.position = "none")
DE_DD_grid_plot <- plot_grid(DE_plot, DD_plot, ncol=1, align="v", labels=c("A", "B"))
multiplot_facet <- ggdraw() +
  draw_plot(DE_DD_grid_plot, 0, 0.1, 1, 0.9) +
  draw_plot(plot_legend, 0, 0, 1, 0.1)
print(multiplot_facet)
dev.off()

# Venn barplot
## DD+ MDSeq non DE genes
all_genes_output_basename <- sprintf("fig5-TCGA_%s_vs_%s_FC_%d_DD_count_barplot", pheno_1, pheno_2, DD_FC_threshold)
df2ggplot <- build_common_specific_barplot_data_frame(all_datasets_MDSeq_genes, all_datasets_DiPhiSeq_DDp_genes, "DDp_non_DE")
DDp_non_DE_barplot <- common_specific_barplot(df2ggplot)
pdf(file=sprintf("%s/%s.pdf", output_dir, all_genes_output_basename), width=9)
print(DDp_non_DE_barplot)
dev.off()

## DD+ MDSeq lowly DE genes
lowly_DE_genes_output_basename <- sprintf("figS7-TCGA_%s_vs_%s_lowly_DE_DD_count_barplot", pheno_1, pheno_2)
df2ggplot <- build_common_specific_barplot_data_frame(all_datasets_MDSeq_genes, all_datasets_DiPhiSeq_DDp_genes, "DDp_lowly_DE")
DDp_lowly_DE_barplot <- common_specific_barplot(df2ggplot)
pdf(file=sprintf("%s/%s.pdf", output_dir, lowly_DE_genes_output_basename), width=9)
print(DDp_lowly_DE_barplot)
dev.off()

