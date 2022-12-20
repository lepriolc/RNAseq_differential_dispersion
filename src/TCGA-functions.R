# DiPhiSeq and MDSeq results reading functions
## convert Ensembl gene IDs to NCBI gene (formerly Entrezgene) IDs
ENSG2entrez <- function(ensg_ids, ensembl, out_dir, out_basename) {
  # remove dot in Ensembl gene IDs
  ensg_no_dot_ids <- unlist(lapply(ensg_ids, function(x) { strsplit(x, "[.]")[[1]][1] }))
  if (length(ensg_no_dot_ids) == 0) {
    entrez_ids <- NA
  } else {
    # get NCBI gene IDs
    gene.data <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = ensg_no_dot_ids, mart = ensembl)
    write.csv(gene.data, file=sprintf("%s/%s.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
    entrez_ids <- unique(gene.data$entrezgene_id[!is.na(gene.data$entrezgene_id)])
  }
  return(entrez_ids)
}

## get DE genes and non-DE genes in MDSeq results
get_MDSeq_DE_genes <- function(mdseq_df, cond_1, cond_2, pval, FC_DE, ensembl, out_dir, out_basename) {
  DE_genes <- list()
  
  ## DE+
  category <- "DE+"
  category_out_dir <- sprintf("%s/%s", out_dir, gsub(" ", "_", category))
  if (! dir.exists(category_out_dir)) {
    dir.create(category_out_dir, mode="0775", recursive=TRUE)
  }
  DEp <- rownames(mdseq_df[which(mdseq_df$FDR.mean < pval & mdseq_df[[sprintf("%svs%s.mean.log2FC.%s", make.names(cond_1), make.names(cond_2), FC_DE)]] > 0),])
  DEp_no_dot <- unlist(lapply(DEp, function(x) { strsplit(x, "[.]")[[1]][1] }))
  DEp_entrez <- ENSG2entrez(DEp_no_dot, ensembl, category_out_dir, sprintf("%s_%s_genes", out_basename, gsub(" ", "_", category)))
  DE_genes[[category]] <- DEp_entrez
  ## DE-
  category <- "DE-"
  category_out_dir <- sprintf("%s/%s", out_dir, gsub(" ", "_", category))
  if (! dir.exists(category_out_dir)) {
    dir.create(category_out_dir, mode="0775", recursive=TRUE)
  }
  DEm <- rownames(mdseq_df[which(mdseq_df$FDR.mean < pval & mdseq_df[[sprintf("%svs%s.mean.log2FC.%s", make.names(cond_1), make.names(cond_2), FC_DE)]] < 0),])
  DEm_no_dot <- unlist(lapply(DEm, function(x) { strsplit(x, "[.]")[[1]][1] }))
  DEm_entrez <- ENSG2entrez(DEm_no_dot, ensembl, category_out_dir, sprintf("%s_%s_genes", out_basename, gsub(" ", "_", category)))
  DE_genes[[category]] <- DEm_entrez
  ## no DE
  category <- "non DE"
  category_out_dir <- sprintf("%s/%s", out_dir, gsub(" ", "_", category))
  if (! dir.exists(category_out_dir)) {
    dir.create(category_out_dir, mode="0775", recursive=TRUE)
  }
  non_DE <- rownames(mdseq_df[which(mdseq_df$FDR.mean > pval),])
  non_DE_no_dot <- unlist(lapply(non_DE, function(x) { strsplit(x, "[.]")[[1]][1] }))
  non_DE_entrez <- ENSG2entrez(non_DE_no_dot, ensembl, category_out_dir, sprintf("%s_%s_genes", out_basename, gsub(" ", "_", category)))
  DE_genes[[category]] <- non_DE_entrez
  
  return(DE_genes)
}

## get DV genes among non-DE genes in MDSeq results
get_MDSeq_DV_no_DE_genes <- function(mdseq_df, cond_1, cond_2, pval, FC_DD, ensembl, out_dir, out_basename) {
  DV_no_DE_genes <- list()
  
  # non-DE DV+
  category <- "non-DE DV+"
  category_out_dir <- sprintf("%s/%s", out_dir, gsub(" ", "_", category))
  if (! dir.exists(category_out_dir)) {
    dir.create(category_out_dir, mode="0775", recursive=TRUE)
  }
  DVp_no_DE <- rownames(mdseq_df[which(mdseq_df$FDR.dispersion < pval & mdseq_df[[sprintf("%svs%s.dispersion.log2FC.%s", make.names(cond_1), make.names(cond_2), FC_DD)]] > 0),])
  DVp_no_DE_no_dot <- unlist(lapply(DVp_no_DE, function(x) { strsplit(x, "[.]")[[1]][1] }))
  DVp_no_DE_entrez <- ENSG2entrez(DVp_no_DE_no_dot, ensembl, category_out_dir, sprintf("%s_%s_genes", out_basename, gsub(" ", "_", category)))
  DV_no_DE_genes[[category]] <- DVp_no_DE_entrez
  
  return(DV_no_DE_genes)
}

## get DD genes in DiPhiseq results
get_DiPhiSeq_DD_genes <- function(diphiseq_df, cond_1, cond_2, pval, ensembl, out_dir, out_basename) {
  pval_colname <- "fdr"
  DD_genes <- list()
  
  ## DD+
  category <- "DD+"
  one_output_dir <- sprintf("%s/%s", out_dir, gsub(" ", "_", category))
  if (! dir.exists(one_output_dir)) {
    dir.create(one_output_dir, mode="0775", recursive=TRUE)
  }
  DDp <- rownames(diphiseq_df[which(diphiseq_df[[sprintf("%s.phi", pval_colname)]] < pval & diphiseq_df[[sprintf("phi.%s", make.names(cond_1))]] > diphiseq_df[[sprintf("phi.%s", make.names(cond_2))]] ),])
  DDp_no_dot <- unlist(lapply(DDp, function(x) { strsplit(x, "[.]")[[1]][1] }))
  DDp_entrez <- ENSG2entrez(DDp_no_dot, ensembl, one_output_dir, sprintf("%s_%s_MDSeq_non_DE_genes", out_basename, gsub(" ", "_", category)))
  DD_genes[[category]] <- DDp_entrez
  
  return(DD_genes)
}

## get DD genes in Levene results
get_Levene_DV_genes <- function(Levene_df, cond_1, cond_2, pval, ensembl, out_dir, out_basename) {
  DV_genes <- list()
  
  # DV+
  category <- "DV+"
  one_output_dir <- sprintf("%s/%s", out_dir, gsub(" ", "_", category))
  if (! dir.exists(one_output_dir)) {
    dir.create(one_output_dir, mode="0775", recursive=TRUE)
  }
  DVp <- rownames(Levene_df[which(Levene_df$p_value.BH < 0.05 & Levene_df[, sprintf("var_%s", cond_1)] > Levene_df[, sprintf("var_%s", cond_2)]),])
  DVp_no_dot <- unlist(lapply(DVp, function(x) { strsplit(x, "[.]")[[1]][1] }))
  DVp_entrez <- ENSG2entrez(DVp_no_dot, ensembl, one_output_dir, sprintf("%s_%s_MDSeq_non_DE_genes", out_basename, gsub(" ", "_", category)))
  DV_genes[[category]] <- DVp_entrez
  
  return(DV_genes)
}

# GO term cluster analysis
## GO term enrichment analysis
go_analysis <- function(gene_list, category, ontology) {
  if (is.null(gene_list)) {
    return(NULL)
  } else {
    return(enrichGO(gene_list[[category]], 'org.Hs.eg.db', ont=ontology))
  }
}

remove_GO_terms_not_in_gosemsim <- function(x, gosemsim) {
  go_terms_in_gosemsim <- names(gosemsim@IC)
  # filter enrichResult object
  if (! is.null(x)) {
    x_df <- as.data.frame(x)
    filt_df <- x_df[rownames(x_df) %in% go_terms_in_gosemsim,]
    x_filt <- x
    slot(x_filt, "result") <- filt_df
  } else {
    x_filt <- x
  }
  
  return(x_filt)
}

similarity_df <- function(enrich, gosemsim, sim_meth) {
  # semantic similarity between GO terms
  GO_IDs <- enrich$ID
  sim_df <- mgoSim(GO_IDs, GO_IDs, semData=gosemsim, measure=sim_meth, combine=NULL)
  return(sim_df)
}

simplify_enrichGO <- function(enrich, gosemsim, sim_meth, sim_thres) {
  if (! is.null(enrich)) {
    # only keep significant enriched GO terms
    enrich_df <- enrich[which(enrich$p.adjust < 0.05),]
    enrich_filt <- enrich
    slot(enrich_filt, "result") <- enrich_df
    # reduce redundancy of enriched GO terms
    filt_simp <- simplify(enrich_filt, cutoff=sim_thres, by="pvalue", measure=sim_meth, semData=gosemsim)
    return(filt_simp$ID)
  } else {
    return(NULL)
  }
}

reduce_redundancy <- function(enrich, gosemsim, sim_meth, sim_thres, output_dir, output_basename) {
  # compute similarity between GO terms for each dataset
  sim_GO <- lapply(enrich, similarity_df, gosemsim, sim_meth)
  ## plots: dendrogram and heatmap
  for (one_dataset in names(sim_GO)) {
    one_sim_df <- sim_GO[[one_dataset]]
    if (! is.null(dim(one_sim_df))) { # one_sim_df is NA (and dim(NA) is null, avoid warning messages for non NA data frames) for datasets without genes in the DE/DV gene category
      nb_GO_terms <- dim(one_sim_df)[1]
      if (nb_GO_terms > 1) { # need at least 2 GO terms to perform clustering
        similarity_output_basename <- sprintf("%s_%s_enrichGO_similarity_%s", one_dataset, output_basename, sim_meth)
        write.csv(one_sim_df, file=sprintf("%s/%s.csv", output_dir, similarity_output_basename), quote=FALSE, row.names=TRUE)
        # cluster GO terms based on 1-similarity
        ## replace NA values by the minimum value (hclust does not handle NA values)
        one_sim_df[is.na(one_sim_df)] <- min(one_sim_df, na.rm=TRUE)
        one_clust <- hclust(as.dist(1-one_sim_df))
        ## plots
        pdf(sprintf("%s/%s.pdf", output_dir, similarity_output_basename))
        ### dendrogram
        if (nb_GO_terms > 2) { # dendogram can only be plotted with at least 3 terms
          plot(one_clust, cex=0.3, main="", xlab="", ylab="1 - similarity", sub="")
          abline(h=1-sim_thres, col="red", lty=2, lwd=2)
        }
        ### heatmap
        pheatmap(one_sim_df, cluster_rows=one_clust, cluster_cols=one_clust)
        dev.off()
        ## GO term list in clusters
        ### get clusters using similarity threshold 
        grp <- cutree(one_clust, h=1-sim_thres)
        ### list GO terms for cluster with more than 1 GO term
        count_per_grp <- table(grp)
        clusters <- as.integer(names(count_per_grp[count_per_grp > 1]))
        if (length(clusters) > 0) {
          original_terms_vector <- c()
          num_cluster_vector <- c()
          for (one_cluster in clusters) {
            go_terms_in_cluster <- names(grp[grp==one_cluster])
            original_terms_vector <- c(original_terms_vector, paste(go_terms_in_cluster, collapse=","))
            num_cluster_vector <- c(num_cluster_vector, one_cluster)
          }
          terms_in_clusters_df <- data.frame(num_cluster_vector, original_terms_vector)
          colnames(terms_in_clusters_df) <- c("num_cluster", "GO_terms")
          write.table(terms_in_clusters_df, file=sprintf("%s/%s_clusters_%s.tsv", output_dir, similarity_output_basename, sub("[.]", "_", sim_thres)), quote=FALSE, sep="\t", row.names=FALSE)
        }
      }
    }
  }
  
  # simplify GO term lists by using semantic similarity between GO terms and a thresold value
  # print("Reduce redundancy among GO terms")
  simp_GO_list <- lapply(enrich, simplify_enrichGO, gosemsim, sim_meth, sim_thres)
  ## write simplified enriched GO term lists
  for (one_dataset in names(simp_GO_list)) {
    one_enrich <- enrich[[one_dataset]]
    if (! is.null(one_enrich)) { # one_enrich is NULL for datasets without genes in the DE/DV gene category
      simp_GO_terms <- simp_GO_list[[one_dataset]]
      write.table(slot(one_enrich, "result")[simp_GO_terms,], file=sprintf("%s/%s_%s_enrichGO_simplified_%s_%s.tsv", output_dir, one_dataset, output_basename, sim_meth, sub("[.]", "_", sim_thres)), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
  
  return(simp_GO_list)
}

# get the common generic term in Gene Ontology of a list of GO terms
GO_common_generic_term <- function(go_terms) {
  nb_terms <- length(go_terms)
  if (nb_terms >= 2) {
    common_generic_term <- getMinimumSubsumer(go_terms[1], go_terms[2])
    if (nb_terms > 2) {
      for (i in 3:nb_terms) {
        if (! common_generic_term %in% getAncestors()$go_terms[i]) {
          common_generic_term <- getMinimumSubsumer(common_generic_term, go_terms[i])
        }
      }
    }
  } else {
    print("Less than 2 GO terms")
    common_generic_term <- NA
  }
  return(common_generic_term)
}

# replace GO terms belonging to a cluster by the corresponding common generic GO term
replace_GO_terms_by_generic_term <- function(x, groups, generic_terms_df) {
  if (x %in% names(groups[groups %in% generic_terms_df$num_cluster])) {
    cluster <- groups[x]
    to_return <- as.character(generic_terms_df[which(generic_terms_df$num_cluster==cluster),"generic_term"])
  } else {
    to_return <- x
  }
  return(to_return)
}

# filter enriched GO term list by extracting a list of GO terms
extract_GO_terms <- function(x, GO_list) {
  if (! is.null(x)) {
    # filter enrichResult object
    x_df <- as.data.frame(x)
    filt_df <- x_df[rownames(x_df) %in% GO_list,]
    x_filt <- x
    slot(x_filt, "result") <- filt_df
  } else {
    x_filt <- x
  }
  
  return(x_filt)
}

GO_cluster_pipeline <- function(enrich, gosemsim, sim_meth, sim_thres, output_dir, output_basename) {
  
  # reduce redundancy in the enriched GO term list of each dataset
  print("Compute similarity between enriched GO terms and reduce redundancy for each dataset")
  reduce_redundancy_output_dir <- sprintf("%s/00-Intra-dataset_redundancy_reduction", output_dir)
  if (! dir.exists(reduce_redundancy_output_dir)) {
    dir.create(reduce_redundancy_output_dir, mode="0775", recursive=TRUE)
  }
  simp_GO_list <- reduce_redundancy(enrich, gosemsim, sim_meth, sim_thres, reduce_redundancy_output_dir, output_basename)
  
  # cluster similar GO terms among all datasets
  print("Cluster similar enriched GO terms among all datasets")
  reduce_redundancy_output_dir <- sprintf("%s/10-Inter-dataset_redundancy_reduction", output_dir)
  if (! dir.exists(reduce_redundancy_output_dir)) {
    dir.create(reduce_redundancy_output_dir, mode="0775", recursive=TRUE)
  }
  GO_list <- unique(unlist(simp_GO_list, use.names=FALSE))
  ## similarity between GO term lists
  one_sim_df <- mgoSim(GO_list, GO_list, semData=gosemsim, measure=sim_meth, combine=NULL)
  similarity_output_basename <- sprintf("TCGA_%s_enrichGO_simpified_%s_%s", output_basename, sim_meth, sub("[.]", "_", sim_thres))
  write.csv(one_sim_df, file=sprintf("%s/%s.csv", reduce_redundancy_output_dir, similarity_output_basename), quote=FALSE, row.names=TRUE)
  ## cluster GO terms based on 1-similarity
  ### replace NA values by the minimum value (hclust does not handle NA values)
  one_sim_df[is.na(one_sim_df)] <- min(one_sim_df, na.rm=TRUE)
  clust <- hclust(as.dist(1-one_sim_df))
  ### plots
  pdf(sprintf("%s/%s.pdf", reduce_redundancy_output_dir, similarity_output_basename))
  #### dendrogram
  if (length(clust$labels) > 2) {
    plot(clust, cex=0.3)
    abline(h=1-sim_thres, col="red", lty=2, lwd=2)
  }
  #### heatmap
  pheatmap(one_sim_df, cluster_rows=clust, cluster_cols=clust)
  dev.off()
  
  # get common generic term for all GO term clusters
  print("Get common generic term for all GO term clusters")
  ## get clusters using similarity threshold 
  grp <- cutree(clust, h=1-sim_thres)
  ## get common generic term of GO terms of clusters with more than 1 GO term
  count_per_grp <- table(grp)
  clusters <- as.integer(names(count_per_grp[count_per_grp > 1]))
  if (length(clusters) > 0) {
    common_generic_terms_vector <- c()
    original_terms_vector <- c()
    num_cluster_vector <- c()
    for (one_cluster in clusters) {
      go_terms_in_cluster <- names(grp[grp==one_cluster])
      ### get common generic GO term
      common_generic_term <- GO_common_generic_term(go_terms_in_cluster)
      ### save terms
      common_generic_terms_vector <- c(common_generic_terms_vector, common_generic_term)
      original_terms_vector <- c(original_terms_vector, paste(go_terms_in_cluster, collapse=","))
      num_cluster_vector <- c(num_cluster_vector, one_cluster)
    }
    common_generic_terms_df <- data.frame(common_generic_terms_vector, original_terms_vector, num_cluster_vector)
    colnames(common_generic_terms_df) <- c("generic_term", "GO_terms", "num_cluster")
    write.table(common_generic_terms_df, file=sprintf("%s/%s_cluster_generic_terms.tsv", reduce_redundancy_output_dir, similarity_output_basename), quote=FALSE, sep="\t", row.names=FALSE)
    ## replace similar GO terms by their common generic term
    GO_list_with_generic_terms <- unlist(lapply(GO_list, replace_GO_terms_by_generic_term, grp, common_generic_terms_df))
  } else {
    GO_list_with_generic_terms <- GO_list
  }
  
  # filter enriched GO term list of each dataset: remove GO terms in similarity clusters and only retain cluster common generic GO terms
  enrichResult_simp_generic_list <- lapply(enrich, extract_GO_terms, GO_list_with_generic_terms)
  for (one_dataset in names(simp_GO_list)) {
    if (! is.null(enrichResult_simp_generic_list[[one_dataset]])) { # enrichResult_simp_generic_list[[one_dataset]] is NULL for datasets without genes of interest
      one_dataset_output_basename <- sprintf("%s_%s", one_dataset, output_basename)
      write.table(slot(enrichResult_simp_generic_list[[one_dataset]], "result"), file=sprintf("%s/%s_%s_enrichGO_simplified_%s_%s_cluster_generic_terms.tsv", reduce_redundancy_output_dir, one_dataset, output_basename, sim_meth, sub("[.]", "_", sim_thres)), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
  
  # results data frame to return
  all_datasets_enrichResult_simp <- merge_result(enrichResult_simp_generic_list)
  df2return <- slot(all_datasets_enrichResult_simp, "compareClusterResult")
  return(df2return)
}

dotplot2customplot <- function(cluster_df, goterms, out_dir, out_basename) {
  # re-order cluster_df by number of datasets with a significant p-value and by p-value mean
  ## order GO IDs by number of datasets with a significant p-value
  datasets_per_GO_ID <- table(cluster_df$ID)
  GO_IDs_order <- names(datasets_per_GO_ID[order(datasets_per_GO_ID, decreasing=TRUE)])
  ## mean -log10(p-value) for each GO ID
  GO_IDs_pval <- unlist(lapply(GO_IDs_order, function(x) { mean(-log10(cluster_df[which(cluster_df$ID==x), "p.adjust"])) }))
  names(GO_IDs_pval) <- GO_IDs_order
  all_datasets <- levels(cluster_df$Cluster)
  GO_IDs_vector <- c()
  mean_pval_vector <- c()
  pval_matrix <- c()
  for (count in unique(datasets_per_GO_ID[order(datasets_per_GO_ID, decreasing=TRUE)])) {
    count_GO_IDs_pval <- GO_IDs_pval[names(datasets_per_GO_ID[which(datasets_per_GO_ID == count)])]
    count_GO_IDs_pval <- count_GO_IDs_pval[order(count_GO_IDs_pval, decreasing=TRUE)]
    ### get p-value per dataset
    for (one_GO_ID in names(count_GO_IDs_pval)) {
      one_sub_df <- cluster_df[which(cluster_df$ID==one_GO_ID),]
      rownames(one_sub_df) <- one_sub_df$Cluster
      pval_matrix <- rbind(pval_matrix, t(one_sub_df[levels(cluster_df$Cluster), "p.adjust"]))
    }
    GO_IDs_vector <- c(GO_IDs_vector, names(count_GO_IDs_pval))
    mean_pval_vector <- c(mean_pval_vector, as.numeric(count_GO_IDs_pval))
  }
  colnames(pval_matrix) <- all_datasets
  rownames(pval_matrix) <- GO_IDs_vector
  
  # data frame for custom plot
  ## number of datasets with significant p-values for each GO ID
  nb_datasets <- datasets_per_GO_ID[GO_IDs_vector]
  ## get description for each GO ID
  description <- goterms[GO_IDs_vector]
  df2customplot <- as.data.frame(cbind(GO_IDs_vector, description, nb_datasets, mean_pval_vector, pval_matrix))
  write.table(df2customplot, file=sprintf("%s/%s_customplot.tsv", out_dir , out_basename), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  
  # customplot
  df2customplot$nb_datasets <- as.numeric(as.character(df2customplot$nb_datasets))
  max_pvalue <- max(-log10(as.numeric(as.character(pval_matrix))), na.rm=TRUE)
  nb_GO_terms_per_page <- 40
  nb_GO_terms_to_plot <- dim(df2customplot)[1]
  i <- 0
  pdf(sprintf("%s/%s_customplot.pdf", out_dir, out_basename), width=9, height=7)
  while(i*nb_GO_terms_per_page < nb_GO_terms_to_plot) {
    ## extract data for custom plot by page
    max_to_plot <- (i+1)*nb_GO_terms_per_page
    if (max_to_plot > nb_GO_terms_to_plot) {
      max_to_plot <- nb_GO_terms_to_plot
    }
    df2customplot_subdf <- df2customplot[(i*nb_GO_terms_per_page+1):max_to_plot,]
    if (dim(df2customplot_subdf)[1] > 0) {
      desc_vector <- c()
      dataset_vector <- c()
      pval_vector <- c()
      for (one_GO in rownames(df2customplot_subdf)) {
        for (one_dataset in all_datasets) {
          desc_vector <- c(desc_vector, as.character(df2customplot_subdf[one_GO, "description"]))
          dataset_vector <- c(dataset_vector, one_dataset)
          pval_vector <- c(pval_vector, -log10(as.numeric(as.character(pval_matrix[one_GO,one_dataset]))))
        }
      }
      df2ggplot <- data.frame(description=desc_vector, dataset=dataset_vector, pvalue=pval_vector)
      df2ggplot$description <- factor(df2ggplot$description, levels=as.character(rev(df2customplot_subdf$description)))
      ## customplot
      p <- ggplot(data = df2ggplot, aes(x=dataset, y=description, fill=pvalue)) + 
        geom_tile() +
        scale_fill_gradient(low="blue", high="red", na.value="white", limits=c(0,max_pvalue)) +
        labs(fill="-log10(p-value)") +
        theme_bw() +
        theme(legend.position = "top") +
        theme(legend.title=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) +
        theme(axis.text.x=element_text(size=12, angle=45, hjust=1), axis.text.y=element_text(size=10, hjust=1), legend.text=element_text(size=9))
      print(p)
    }
    i <- i+1
  }
  dev.off()
}
