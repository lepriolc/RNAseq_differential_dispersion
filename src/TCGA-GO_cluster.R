library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(pheatmap)
library(GO.db)
library(GOSim)
library(ggplot2)
library(gridExtra)
library(biomaRt)


##############
# Parameters #
##############

# datasets
datasets <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-THCA")
# conditions
pheno_1 <- "TP"
pheno_2 <- "NT"
# fold-change threshold to identify DE genes (not in log2 scale)
DE_FC_threshold <- 1
# fold-change threshold to identify DD genes (not in log2 scale)
DD_FC_threshold <- 1
# p-value threshold
pval_threshold <- 0.05
# normalization method
normalization_method <- "TMM"
# GO term enrichment analysis
## ontology
ontology <- "BP"
## similarity method
similarity_method <- "Rel"
## similarity threshold
similarity_threshold <- 0.8
# path to output directory
output_dir <- "./output/TCGA"


#############
# Functions #
#############
source("./src/TCGA-functions.R")


#############
# Load data #
#############
# select Human genes (GRCh38.p13) dataset in BioMart ensemsbl database
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# get all GO terms/IDs
goterms <- Term(GOTERM)


############
# Analysis #
############
# get lowly DE genes and DD genes among lowly DE genes
go_cluster_output_dir <- sprintf("%s/60-GO_cluster", output_dir)
gene_list_output_dir <- sprintf("%s/00-Gene_lists", go_cluster_output_dir) # TODO before publication: replace by next line

all_datasets_MDSeq_DE_genes <- list()
all_datasets_DV_DD_non_DE_genes <- list()

## read Levene's test, MDSeq, DiPhiSeq, GAMLSS and DiffDist output files
for (one_dataset in datasets) {
  print(one_dataset)
  # DE genes according to MDSeq
  mdseq_output_dir <- sprintf("%s/20-MDSeq/%s", output_dir, one_dataset)
  mdseq_output_basename <- sprintf("%s_%s_vs_%s_%s_MDSeq_results", one_dataset, pheno_1, pheno_2, normalization_method)
  mdseq_de_file <- sprintf("%s/%s.csv", mdseq_output_dir, mdseq_output_basename)
  mdseq_de_df <- read.csv(file=mdseq_de_file, header=TRUE, row.names=1)
  ## get Entrez IDs
  mdseq_gene_list_output_dir <- sprintf("%s/20-MDSeq", gene_list_output_dir)
  all_datasets_MDSeq_DE_genes[[one_dataset]] <- get_MDSeq_DE_genes(mdseq_de_df, pheno_1, pheno_2, pval_threshold, log(DE_FC_threshold, 2), ensembl, mdseq_gene_list_output_dir, mdseq_output_basename)
  ## get non-DE genes
  mdseq_non_DE_genes <- all_datasets_MDSeq_DE_genes[[one_dataset]][["non DE"]]
  
  # DD genes in MDSeq non-DE genes
  ## MDSeq
  mdseq_dv_df <- mdseq_de_df[which(mdseq_de_df$FDR.mean > pval_threshold),]
  all_datasets_DV_DD_non_DE_genes[[one_dataset]] <- get_MDSeq_DV_no_DE_genes(mdseq_dv_df, pheno_1, pheno_2, pval_threshold, log(DD_FC_threshold, 2), ensembl, mdseq_gene_list_output_dir, mdseq_output_basename)
  
  ## Levene's test
  Levene_output_dir <- sprintf("%s/10-Levene/%s", output_dir, one_dataset)
  Levene_output_basename <- sprintf("%s_%s_vs_%s_%s_Levene_test_results", one_dataset, pheno_1, pheno_2, normalization_method)
  Levene_results_file <- sprintf("%s/%s.csv", Levene_output_dir, Levene_output_basename)
  Levene_results_df <- read.csv(file=Levene_results_file, header=TRUE, row.names=1)
  ### get DD genes Entrez IDs
  Levene_gene_list_output_dir <- sprintf("%s/10-Levene", gene_list_output_dir)
  Levene_DV_genes <- get_Levene_DV_genes(Levene_results_df, pheno_1, pheno_2, pval_threshold, ensembl, Levene_gene_list_output_dir, Levene_output_basename)
  ### get DD+ genes in MDSeq non-DE genes
  mdseq_category <-"non-DE DV+"
  Levene_DVp_genes <- Levene_DV_genes[["DV+"]]
  Levene_DVp_MDSeq_non_DE_genes <- Levene_DVp_genes[Levene_DVp_genes %in% mdseq_non_DE_genes]
  Levene_DVp_MDSeq_non_DE_genes_specific <- Levene_DVp_MDSeq_non_DE_genes[! Levene_DVp_MDSeq_non_DE_genes %in% all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]]]
  ### add Levene specific DD+ genes in MDSeq no DE genes
  all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]] <- c(all_datasets_DV_DD_non_DE_genes[[one_dataset]][[mdseq_category]], Levene_DVp_MDSeq_non_DE_genes_specific)
  
  ## DiPhiSeq
  diphiseq_output_dir <- sprintf("%s/30-DiPhiSeq/%s", output_dir, one_dataset)
  diphiseq_output_basename <- sprintf("%s_%s_vs_%s_%s_DiPhiSeq_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diphiseq_results_file <- sprintf("%s/%s.csv", diphiseq_output_dir, diphiseq_output_basename)
  diphiseq_results_df <- read.csv(file=diphiseq_results_file, header=TRUE, row.names=1)
  ### get DD genes Entrez IDs
  diphiseq_gene_list_output_dir <- sprintf("%s/30-DiPhiSeq", gene_list_output_dir)
  DiPhiSeq_DD_genes <- get_DiPhiSeq_DD_genes(diphiseq_results_df, pheno_1, pheno_2, pval_threshold, ensembl, diphiseq_gene_list_output_dir, diphiseq_output_basename)
  ### get DD+ genes in MDSeq non-DE genes
  mdseq_category <-"non-DE DV+"
  diphiseq_DDp_genes <- DiPhiSeq_DD_genes[["DD+"]]
  diphiseq_DDp_MDSeq_non_DE_genes <- diphiseq_DDp_genes[diphiseq_DDp_genes %in% mdseq_non_DE_genes]
  diphiseq_DDp_MDSeq_non_DE_genes_specific <- diphiseq_DDp_MDSeq_non_DE_genes[! diphiseq_DDp_MDSeq_non_DE_genes %in% all_datasets_DV_DD_non_DE_genes[[one_dataset]][[mdseq_category]]]
  ### add DiPhiSeq specific DD+ genes in MDSeq no DE genes
  all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]] <- c(all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]], diphiseq_DDp_MDSeq_non_DE_genes_specific)
  
  ## GAMLSS
  gamlss_output_dir <- sprintf("%s/40-GAMLSS/%s", output_dir, one_dataset)
  gamlss_output_basename <- sprintf("%s_%s_vs_%s_%s_GAMLSS_results", one_dataset, pheno_1, pheno_2, normalization_method)
  gamlss_results_file <- sprintf("%s/%s.csv", gamlss_output_dir, gamlss_output_basename)
  gamlss_results_df <- read.csv(file=gamlss_results_file, header=TRUE, row.names=1)
  ### get DD genes
  gamlss_DDp_genes <- rownames(gamlss_results_df[which(gamlss_results_df$padj.cv < pval_threshold & gamlss_results_df[, sprintf("CV.%s", pheno_1)] > gamlss_results_df[, sprintf("CV.%s", pheno_2)]),])
  #### only keep genes whose logFC is consistent with DiPhiSeq log2FC
  gamlss_DD_genes_logFC <- log(gamlss_results_df[, sprintf("CV.%s", pheno_1)]/gamlss_results_df[, sprintf("CV.%s", pheno_2)], 2)
  gamlss_DD_genes_DiPhiSeq_log2FC <- diphiseq_results_df[gamlss_DDp_genes, sprintf("phi.%s", pheno_1)] - diphiseq_results_df[gamlss_DDp_genes, sprintf("phi.%s", pheno_2)]
  gamlss_DDp_genes_consistent_logFC <- c()
  for (i in 1:length(gamlss_DDp_genes)) {
    if (!is.na(gamlss_DD_genes_DiPhiSeq_log2FC[i]) & (gamlss_DD_genes_logFC[i] > 0 & gamlss_DD_genes_DiPhiSeq_log2FC[i] > 0)) {
      gamlss_DDp_genes_consistent_logFC <- c(gamlss_DDp_genes_consistent_logFC, gamlss_DDp_genes[i])
    }
  }
  print(sprintf("%d out of %d GAMLSS DDp genes with consistent FC with DiPhiSeq", length(gamlss_DDp_genes_consistent_logFC), length(gamlss_DDp_genes)))
  gamlss_DDp_genes <- gamlss_DDp_genes_consistent_logFC
  ### get DD genes Entrez IDs
  gamlss_gene_list_output_dir <- sprintf("%s/40-GAMLSS", gene_list_output_dir)
  category <- "DD+"
  category_output_dir <- sprintf("%s/%s", gamlss_gene_list_output_dir, gsub(" ", "_", category))
  if (! dir.exists(category_output_dir)) {
    dir.create(category_output_dir, mode="0775", recursive=TRUE)
  }
  gamlss_DDp_genes_no_dot <- unlist(lapply(gamlss_DDp_genes, function(x) { strsplit(x, "[.]")[[1]][1] }))
  gamlss_DDp_genes_entrez <- ENSG2entrez(gamlss_DDp_genes_no_dot, ensembl, category_output_dir, sprintf("%s_%s_MDSeq_non_DE_genes", gamlss_output_basename, gsub(" ", "_", category)))
  ### get DD+ genes in MDSeq non-DE genes
  gamlss_DDp_MDSeq_non_DE_genes <- gamlss_DDp_genes_entrez[gamlss_DDp_genes_entrez %in% mdseq_non_DE_genes]
  gamlss_DDp_MDSeq_non_DE_genes_specific <- gamlss_DDp_MDSeq_non_DE_genes[! gamlss_DDp_MDSeq_non_DE_genes %in% all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]]]
  ### add GAMLSS specific DD+ genes in MDSeq no DE genes
  all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]] <- c(all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]], gamlss_DDp_MDSeq_non_DE_genes_specific)
  
  ## DiffDist
  diffdist_output_dir <- sprintf("%s/50-DiffDist/%s", output_dir, one_dataset)
  diffdist_output_basename <- sprintf("%s_%s_vs_%s_%s_DiffDist_results", one_dataset, pheno_1, pheno_2, normalization_method)
  diffdist_results_file <- sprintf("%s/%s.csv", diffdist_output_dir, diffdist_output_basename)
  diffdist_results_df <- read.csv(file=diffdist_results_file, header=TRUE, row.names=1)
  ### get DD genes
  diffdist_DDp_genes <- rownames(diffdist_results_df[which(diffdist_results_df[[sprintf("disp.%svs%s.pval.BH", pheno_1, pheno_2)]] < pval_threshold & diffdist_results_df[[sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]] > 0),])
  #### only keep genes whose logFC is consistent with DiPhiSeq log2FC
  diffdist_DD_genes_logFC <- diffdist_results_df[diffdist_DDp_genes, sprintf("disp.%svs%s.logFC", pheno_1, pheno_2)]
  diffdist_DD_genes_DiPhiSeq_log2FC <- diphiseq_results_df[diffdist_DDp_genes, sprintf("phi.%s", pheno_1)] - diphiseq_results_df[diffdist_DDp_genes, sprintf("phi.%s", pheno_2)]
  diffdist_DDp_genes_consistent_logFC <- c()
  for (i in 1:length(diffdist_DDp_genes)) {
    if (!is.na(diffdist_DD_genes_DiPhiSeq_log2FC[i]) & (diffdist_DD_genes_logFC[i] > 0 & diffdist_DD_genes_DiPhiSeq_log2FC[i] > 0)) {
      diffdist_DDp_genes_consistent_logFC <- c(diffdist_DDp_genes_consistent_logFC, diffdist_DDp_genes[i])
    }
  }
  print(sprintf("%d out of %d DiffDist DDp genes with consistent FC with DiPhiSeq", length(diffdist_DDp_genes_consistent_logFC), length(diffdist_DDp_genes)))
  diffdist_DDp_genes <- diffdist_DDp_genes_consistent_logFC
  ### get DD genes Entrez IDs
  diffdist_gene_list_output_dir <- sprintf("%s/50-DiffDist", gene_list_output_dir)
  category <- "DD+"
  category_output_dir <- sprintf("%s/%s", diffdist_gene_list_output_dir, gsub(" ", "_", category))
  if (! dir.exists(category_output_dir)) {
    dir.create(category_output_dir, mode="0775", recursive=TRUE)
  }
  diffdist_DDp_genes_no_dot <- unlist(lapply(diffdist_DDp_genes, function(x) { strsplit(x, "[.]")[[1]][1] }))
  diffdist_DDp_genes_entrez <- ENSG2entrez(diffdist_DDp_genes_no_dot, ensembl, category_output_dir, sprintf("%s_%s_MDSeq_non_DE_genes", diffdist_output_basename, gsub(" ", "_", category)))
  ### get DD+ genes in MDSeq non-DE genes
  diffdist_DDp_MDSeq_non_DE_genes <- diffdist_DDp_genes_entrez[diffdist_DDp_genes_entrez %in% mdseq_non_DE_genes]
  diffdist_DDp_MDSeq_non_DE_genes_specific <- diffdist_DDp_MDSeq_non_DE_genes[! diffdist_DDp_MDSeq_non_DE_genes %in% all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]]]
  ### add DiffDist specific DD+ genes in MDSeq no DE genes
  all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]] <- c(all_datasets_DV_DD_non_DE_genes[[one_dataset]][["non-DE DD+"]], diffdist_DDp_MDSeq_non_DE_genes_specific)
}

# GO term cluster analysis
categories <- c("non-DE DD+", "DE+", "DE-")
# categories <- c("non-DE DD+")
for (one_category in categories) {
  print(sprintf("category: %s", one_category))
  category_output_dir <- sprintf("%s/10-GO_analysis/%s", go_cluster_output_dir, gsub(" ", "_", one_category))
  if (! dir.exists(category_output_dir)) {
    dir.create(category_output_dir, mode="0775", recursive=TRUE)
  }
  
  ## get gene list
  if (one_category %in% c("DE+", "DE-")) {
    gene_list <- all_datasets_MDSeq_DE_genes
  } else {
    gene_list <- all_datasets_DV_DD_non_DE_genes
  }
  
  ## GO term enrichment analysis
  print("GO term enrichment")
  go_term_enrichment_output_dir <- sprintf("%s/00-Enrichment", category_output_dir)
  if (! dir.exists(go_term_enrichment_output_dir)) {
    dir.create(go_term_enrichment_output_dir, mode="0775", recursive=TRUE)
  }
  output_basename <- sprintf("%s_vs_%s_%s", pheno_1, pheno_2, normalization_method)
  enrichResult_list <- lapply(gene_list, go_analysis, one_category, ontology)
  ### write results of GO term enrichment analysis
  for (one_dataset in names(enrichResult_list)) {
    one_enrich <- enrichResult_list[[one_dataset]]
    if (! is.null(one_enrich)) { # one_enrich is NULL for datasets without genes in the DE/DV gene category
      go_term_enrichment_output_basename <- sprintf("%s_%s", one_dataset, output_basename)
      write.table(slot(one_enrich, "result"), file=sprintf("%s/%s_enrichGO.tsv", go_term_enrichment_output_dir, go_term_enrichment_output_basename), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
  
  if (length(enrichResult_list) > 1) {
    ## remove GO terms which are not valid for Human
    hsGOSemSim_withIC <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
    ### remove GO terms which are not in GOSemSimDATA object (i.e. GO terms which are not valid for Human: https://support.bioconductor.org/p/115680/)
    enrichResult_list_only_human <- lapply(enrichResult_list, remove_GO_terms_not_in_gosemsim, hsGOSemSim_withIC)
    ### write enriched human GO terms
    for (one_dataset in names(enrichResult_list_only_human)) {
      one_enrich <- enrichResult_list_only_human[[one_dataset]]
      if (! is.null(one_enrich)) { # one_enrich is NULL for datasets without genes in the DE/DV gene category
        go_term_enrichment_output_basename <- sprintf("%s_%s", one_dataset, output_basename)
        write.table(slot(one_enrich, "result"), file=sprintf("%s/%s_enrichGO_human_enriched.tsv", go_term_enrichment_output_dir, go_term_enrichment_output_basename), sep="\t", quote=FALSE, row.names=FALSE)
      }
    }
    
    ## cluster enriched GO terms by similarity
    cluster_output_dir <- sprintf("%s/10-Cluster", category_output_dir)
    if (! dir.exists(cluster_output_dir)) {
      dir.create(cluster_output_dir, mode="0775", recursive=TRUE)
    }
    print(sprintf("similarity method: %s, similarity threshold: %.1f", similarity_method, similarity_threshold))
    clusterProfiler_df <- GO_cluster_pipeline(enrichResult_list_only_human, hsGOSemSim_withIC, similarity_method, similarity_threshold, cluster_output_dir, output_basename)
    go_cluster_pipeline_output_basename <- sprintf("TCGA_%s_enrichGO_simplified_%s_%s_cluster_generic_terms", output_basename, similarity_method, sub("[.]", "_", similarity_threshold))
    write.table(clusterProfiler_df, file=sprintf("%s/%s.tsv", cluster_output_dir, go_cluster_pipeline_output_basename), sep="\t", quote=FALSE, row.names=FALSE)
    ### reformat dotplot
    dotplot2customplot(clusterProfiler_df, goterms, cluster_output_dir, go_cluster_pipeline_output_basename)
  }
}

