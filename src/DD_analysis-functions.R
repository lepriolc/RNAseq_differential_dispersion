# DD gene identification
filter_lowly_expressed_genes <- function(count_matrix, con, norm, threshold) {
  print(sprintf("Filter lowly expressed genes with CPM threshold: %s", threshold))
  y <- DGEList(counts=count_matrix, group=as.factor(con))
  y <- calcNormFactors(y, method=norm)
  keep_genes <- rowMeans(cpm(y)) >= threshold
  filtered_count_matrix <- count_matrix[keep_genes,]
  return(filtered_count_matrix)
}

run_Levene <- function(count_matrix, conditions, samples_1, samples_2, norm) {
  # order samples
  count_matrix <- count_matrix[, c(samples_1, samples_2)]
  
  # TMM normalization
  y <- DGEList(counts=count_matrix, group=conditions)
  y <- calcNormFactors(y, method=norm)
  ## mean of the normalized library size
  mean_normalized_library_size <- mean(y$samples$lib.size * y$samples$norm.factors)
  norm_matrix <- sapply(1:ncol(count_matrix),function(x) { (count_matrix[,x] / (y$samples[x,"lib.size"]*y$samples[x,"norm.factors"])) * mean_normalized_library_size })
  colnames(norm_matrix) <- colnames(count_matrix)
  rownames(norm_matrix) <- rownames(count_matrix)
  
  # log-transformation
  log_norm_matrix <- log(norm_matrix + 1, 2)
  
  # Levene test
  results <- as.data.frame(t(apply(log_norm_matrix, 1, function(x, groups=conditions, group_samples_1=samples_1, group_samples_2=samples_2) { 
    var1 <- var(x[group_samples_1])
    var2 <- var(x[group_samples_2])
    test_output <- levene.test(x, groups, location = "mean")
    return(c(var_1=var1, var_2=var2, statistic=unname(test_output$statistic), p_value=test_output$p.value))
  })))
  results <- cbind(results, p_value.BH=p.adjust(results$p_value, "BH"), p_value.BY=p.adjust(results$p_value, "BY"))
  colnames(results)[which(colnames(results) == "var_1")] <- sprintf("var_%s", levels(conditions)[1])
  colnames(results)[which(colnames(results) == "var_2")] <- sprintf("var_%s", levels(conditions)[2])
  
  return(results)
}

run_MDSeq_fit <- function(count_matrix, con, outlier_bool, X_cov, U_cov, cpus, min_size, out_dir, out_basename) {
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

run_MDSeq_DD <- function(fit, con, FC) {
  # testing with a given log2 fold change
  one_FC <- log(FC, 2)
  print(sprintf("Get DE and DV genes with log2FC=%.2f", one_FC))
  results <- extract.ZIMD(fit, compare=list(A=levels(con)[1], B=levels(con)[2]), log2FC.threshold=one_FC)
  ## add BH FDR
  FDR.BH.mean <- p.adjust(results$Pvalue.mean, "BH")
  results <- cbind(results, FDR.BH.mean)
  FDR.BH.dispersion <- p.adjust(results$Pvalue.dispersion, "BH")
  results <- cbind(results, FDR.BH.dispersion)
  return(results)
}

run_DiPhiSeq <- function(count_matrix, con, norm, threshold) { # unused threshold parameter
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

run_GAMLSS <- function(count_matrix, conditions, condition_ref, norm) {
  # normalization
  ## edgeR-style code to get lowly expressed genes to filter after normalize read counts and filter lowly expressed genes (discrepancies between edgeR and MDSeq filter functions)
  y <- DGEList(counts=count_matrix, group=as.factor(conditions))
  if (norm == "TMM") {
    ### TMM normalization
    print("TMM normalization")
    y <- calcNormFactors(y, method=norm)
  }
  ## Calculate counts per million (CPM).
  y$CPM <- cpm.DGEList(y)
  ## Calculate offset variable for each library as natural log of the product of library size and normalization factor.
  ofs <- log(y$samples$lib.size*y$samples$norm.factors)
  y$samples$offset <- ofs
  
  # hereafter, code modified from https://github.com/Vityay/ExpVarQuant/blob/master/ExpVarQuant.R
  
  # Uncomment this code when working in Rstudio to prevent it from locking up
  # sink(file = "diversionOfOutput.txt", append = TRUE)
  
  gene_i <- seq_along(y$counts[,1])
  
  # To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) to estimate GAMLSS models for the first hundred genes.
  results <- lapply(gene_i, function(i) {
    
    # For each gene (i) a table containing: x - a factor of interest (age); y - RNA-seq. counts and offset (ofs) is created.
    dat <- data.frame(
      x = y$samples$group,
      y = y$counts[i,],
      ofs = y$samples$offset
    )
    # indicate reference level
    dat$x <- relevel(dat$x, ref = c(condition_ref))
    
    # Fit negative binomial (family = NBI()) GAMLSS model, which accounts for age effects on mean and overdispersion (non-Poisson noise).
    # fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
    # In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
    m0 <- tryCatch(
      gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
             family = NBI(), sigma.start = 0.1, n.cyc = 100),
      warning= function(w) NULL, error= function(e) NULL
    )
    
    # Fit reduced model by omitting age factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
    m1 <- tryCatch(
      gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
             family = NBI(), sigma.start = 0.1, n.cyc = 100),
      warning= function(w) NULL, error= function(e) NULL
    )
    
    # Create data frame results to store the results.
    results <- data.frame(
      cpm.1 = NA,
      cpm.2=NA,
      LR.cpm = NA,
      # p_gamlss.cpm = NA,
      # p_glm.cpm = NA,
      CV.1 = NA,
      CV.2 = NA,
      LR.cv = NA,
      p.cv = NA
    )
    
    # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
    if(!any(sapply(list(m0,m1), is.null))) 
    {
      # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
      results$cpm.1 = exp(m0$mu.coefficients+log(1e06))[[1]]
      results$cpm.2 = exp(m0$mu.coefficients+log(1e06))[[2]]
      
      # Calculate log2 ratio for changes in CPMs between old and young mice.
      results$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])
      
      # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
      results$CV.1  = sqrt(exp(m0$sigma.coefficients)[[1]])
      results$CV.2 = sqrt(exp(m0$sigma.coefficients)[[2]])
      
      # Calculate log2 ratio for changes in cv(μ) between old and young mice.
      results$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))
      
      # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
      # p.cv – p value of LR test statistic: D_α=-2loga[L(μ_j,α_0  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
      results$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
    }
    results
  })
  # closeAllConnections()
  
  # Transform list results containing GAMLSS estimations for each gene to data frame
  results <- do.call(rbind, results)
  rownames(results) <- rownames(y$counts)[gene_i]
  
  # Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
  # First, remove genes, for which GAMLSS model has failed.
  results <- na.omit(results)
  
  # # Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
  # idx <- results$CV.young > 3 | results$CV.young < 1e-03 | results$CV.old > 3 | results$CV.old < 1e-03
  # results <- results[!idx,]
  
  # Finally, calculate false discovery rates to account for multiple hypothesis testing.
  results$padj.cv <- p.adjust(results$p.cv, "fdr")
  results$padj.BY.cv <- p.adjust(results$p.cv, "BY")
  
  ### add group names in column names
  colnames(results)[which(grepl("1", colnames(results)))] <- sub("1", condition_ref, colnames(results)[which(grepl("1", colnames(results)))])
  colnames(results)[which(grepl("2", colnames(results)))] <- sub("2", levels(conditions)[levels(conditions) != condition_ref], colnames(results)[which(grepl("2", colnames(results)))])
  
  return(results)
}

run_DiffDist <- function(count_matrix, conditions, norm) {
  # code modified from https://github.com/aedanr/DiffDist/blob/main/methods_comparisons.R
  
  # normalization
  libsizes <- colSums(count_matrix)
  nf <- calcNormFactors(count_matrix, method=norm)
  els <- nf * libsizes
  sf <- els / exp(mean(log(libsizes)))
  norm.counts <- t(t(count_matrix) / sf)
  
  # hierarchical model
  HM <- ln_hmm_adapt_3_chains(counts=norm.counts, groups=conditions)
  ## differential expression
  mean.diff.HM <- unname(log(as.matrix(HM$means1)) - log(as.matrix(HM$means2)))
  ### posterior samples of between-group differences in means
  p.mean.HM <- apply(mean.diff.HM, 2, hpd.pval) # tail probabilities of no differences in means
  ### mean log fold-change
  mean.logFC.HM <- apply(log(as.matrix(HM$means1)) - log(as.matrix(HM$means2)), 2, mean)
  ## differential dispersion
  disp.diff.HM <- unname(log(as.matrix(HM$disps1)) - log(as.matrix(HM$disps2)))
  ### posterior samples of between-group differences in dispersions
  p.disp.HM <- apply(disp.diff.HM, 2, hpd.pval) # tail probabilities of no differences in dispersions
  ### dispersion log fold-change
  disp.logFC.HM <- apply(log(as.matrix(HM$disps1)) - log(as.matrix(HM$disps2)), 2, mean)
  
  # result data frame
  results <- data.frame(mean.logFC=mean.logFC.HM, mean.pval=p.mean.HM, mean.pval.BH=p.adjust(p.mean.HM, "BH"), mean.pval.BY=p.adjust(p.mean.HM, "BY"),
                        disp.logFC=disp.logFC.HM, disp.pval=p.disp.HM, disp.pval.BH=p.adjust(p.disp.HM, "BH"), disp.pval.BY=p.adjust(p.disp.HM, "BY"))
  rownames(results) <- rownames(count_matrix)
  ## add group names in column names
  colnames(results)[which(grepl("mean", colnames(results)))] <- sub("mean", sprintf("mean.%svs%s", levels(conditions)[1], levels(conditions)[2]), colnames(results)[which(grepl("mean", colnames(results)))])
  colnames(results)[which(grepl("disp", colnames(results)))] <- sub("disp", sprintf("disp.%svs%s", levels(conditions)[1], levels(conditions)[2]), colnames(results)[which(grepl("disp", colnames(results)))])
  
  return(results)
}
