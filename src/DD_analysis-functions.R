# DD gene identification
filter_lowly_expressed_genes <- function(count_matrix, con, norm, threshold) {
  print(sprintf("Filter lowly expressed genes with CPM threshold: %s", threshold))
  y <- DGEList(counts=count_matrix, group=as.factor(con))
  y <- calcNormFactors(y, method=norm)
  keep_genes <- rowMeans(cpm(y)) >= threshold
  filtered_count_matrix <- count_matrix[keep_genes,]
  return(filtered_count_matrix)
}

run_DiPhiSeq <- function(count_matrix, con, norm, threshold) {
  # normalization
  ## get estimated sequencing depths
  y <- DGEList(counts=count_matrix, group=as.factor(con))
  if (norm == "TMM") {
    ## TMM normalization
    print("TMM normalization")
    y <- calcNormFactors(y, method=norm)
  }
  est_seq_depth <- y$samples$lib.size * y$samples$norm.factors
  
  # DiPhiSeq
  print("DiPhiSeq")
  results <- diphiseq(count_matrix, con, depth=est_seq_depth)
  
  return(list(results=results, samples=y$samples))
}

MDSeq_fit <- function(count_matrix, con, outlier_bool, X_cov, U_cov, cpus, min_size, out_dir, out_basename) {
  design <- get.model.matrix(con)
  if (outlier_bool) {
    # remove outliers
    print("Remove outliers for mean and dispersion")
    count.checked <- remove.outlier(count_matrix, contrast=list(mean=design$mean, dispersion=design$dispersion), mc.cores=cpus, min.sample.size=min_size)
    write.csv(count.checked$outliers, file=sprintf("%s/%s_outliers.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
    counts <- count.checked$count[count.checked$outliers$status==0,]
  } else {
    counts <- count_matrix
  }
  # fit GLM for DE and DV genes
  print("Fit GLM for DE and DV genes")
  fit <- MDSeq(counts, X=X_cov, U=U_cov, contrast=list(mean=design$mean, dispersion=design$dispersion), mc.cores=cpus)
  return(fit)
}

MDSeq_DE_DD <- function(fit, con, FC) {
  # testing with a given log2 fold change
  one_FC <- log(FC, 2)
  print(sprintf("Get DE and DV genes with log2FC=%.2f", one_FC))
  mdseq_results <- extract.ZIMD(fit, compare=list(A=levels(con)[1], B=levels(con)[2]), log2FC.threshold=one_FC)
  ## add BH FDR
  FDR.BH.mean <- p.adjust(mdseq_results$Pvalue.mean, "BH")
  mdseq_results <- cbind(mdseq_results, FDR.BH.mean)
  FDR.BH.dispersion <- p.adjust(mdseq_results$Pvalue.dispersion, "BH")
  mdseq_results <- cbind(mdseq_results, FDR.BH.dispersion)
  return(mdseq_results)
}

MDSeq_DD_for_lowly_DE <- function(count_matrix, con, norm, outlier_bool, X_cov, U_cov, cores, min_size, FC_DE, FC_DD, pval, out_dir, out_basename) {
  # normalization
  exp.normalized <- normalize.counts(count_matrix, group=con, method=norm)
  # fit MDSeq model
  out_basename <- sprintf("%s_%s_MDSeq", out_basename, norm)
  fit <- MDSeq_fit(exp.normalized, con, outlier_bool, X_cov, U_cov, cores, min_size, out_dir, out_basename)
  write.csv(fit$Dat, file=sprintf("%s/%s_fit_Dat.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
  if (! is.null(X_cov)) {
    write.csv(fit$X, file=sprintf("%s/%s_fit_X.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
  }
  if (! is.null(U_cov)) {
    write.csv(fit$U, file=sprintf("%s/%s_fit_U.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
  }
  # run MDSeq to filer highly DE genes
  mdseq_DE_results <- MDSeq_DE_DD(fit, con, FC_DE)
  write.csv(mdseq_DE_results, file=sprintf("%s/%s_FC_%s.csv", out_dir, out_basename, sub("[.]", "_", FC_DE)), quote=FALSE, row.names=TRUE)
  ## get lowly-DE genes
  lowly_DE_genes <- rownames(mdseq_DE_results[which(mdseq_DE_results$FDR.mean > pval),])
  lowly_DE_counts <- fit$counts[lowly_DE_genes,]
  # fit MDSeq model 
  out_basename <- sprintf("%s_lowly_DE", out_basename)
  fit_lowly_DE <- MDSeq_fit(lowly_DE_counts, con, FALSE, X_cov, U_cov, cores, min_size, out_dir, out_basename)
  write.csv(fit_lowly_DE$Dat, file=sprintf("%s/%s_fit_Dat.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
  if (! is.null(X_cov)) {
    write.csv(fit_lowly_DE$X, file=sprintf("%s/%s_fit_X.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
  }
  if (! is.null(U_cov)) {
    write.csv(fit_lowly_DE$U, file=sprintf("%s/%s_fit_U.csv", out_dir, out_basename), quote=FALSE, row.names=TRUE)
  }
  # run MDSeq to identify DD genes among lowly DE genes
  mdseq_DD_results <- MDSeq_DE_DD(fit_lowly_DE, con, FC_DD)
  return(mdseq_DD_results)
}
