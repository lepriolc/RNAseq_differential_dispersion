# generate datasets
create_effect_size_vector <- function(length, pct, pct_down, base, param, runif_bool) {
  # positive genes (i.e. DE or DD genes) effect
  nb_pos_genes <- floor(pct * length)
  if (param > 0) {
    ## exponential distribution
    effect.size <- base + rexp(nb_pos_genes, param)
  } else {
    ## a constant
    effect.size <- rep(base, nb_pos_genes)
  }
  ## inverse for the values for the pct_down genes
  genes_down <- floor(pct_down * nb_pos_genes)
  if (genes_down > 0) {
    effect.size[1:genes_down] <- unlist(lapply(effect.size[1:genes_down], function(x) { 1/x }))
  }
  
  # negative genes (i.e. non-DE or non-DD genes) effect
  nb_neg_genes <- length - nb_pos_genes
  if (runif_bool) {
    # effect size = runif() with min=1 and max=base effect
    neg.effect.size <- runif(nb_neg_genes, 1, base)
    ## inverse for the values for half of the effect sizes
    # genes_down <- floor(nb_neg_genes/2)
    # neg.effect.size[1:genes_down] <- unlist(lapply(neg.effect.size[1:genes_down], function(x) { 1/x }))
    genes_down <- floor(pct_down * nb_neg_genes)
    if (genes_down > 0) {
      neg.effect.size[1:genes_down] <- unlist(lapply(neg.effect.size[1:genes_down], function(x) { 1/x }))
    }
    effect.size <- c(effect.size, neg.effect.size)
  } else {
    # effect size = 1 for negative genes
    effect.size <- c(effect.size, rep(1, nb_neg_genes))
  }
  
  return(effect.size)
}

get_DE_effect_param <- function(DE_base) {
  if (DE_base == 1.1) {
    DE_param <- 0.85
  } else {
    if (DE_base == 1.2 | DE_base == 1.3) {
      DE_param <- 0.9
    } else {
      if (DE_base == 1.4) {
        DE_param <- 0.95
      } else {
        if (DE_base == 1.5) {
          DE_param <- 1
        }
      }
    }
  }
  return(DE_param)
}

# DD analysis performance
DD_category <- function(FC, FC_threshold, to_print) {
  if (FC > FC_threshold) {
    category <- sprintf("%s+", to_print)
  } else {
    if (FC < -FC_threshold) {
      category <- sprintf("%s-", to_print)
    } else {
      category <- sprintf("non%s", to_print)
    }
  }
  return(category)
}

get_contingency_category <- function(pval, label, threshold) {
  if (is.na(pval)) {
    category <- NA
  } else {
    if (label == 1) {
      # positive
      if (pval < threshold) {
        category <- "TP"
      } else {
        category <- "FN"
      }
    } else {
      if (label == 0) {
        # negative
        if (pval < threshold) {
          category <- "FP"
        } else {
          category <- "TN"
        }
      }
    }
  }
  return(category)
}

add_contingency_category <- function(results_df, pval_colname, label_colname, threshold) {
  category_vector <- c()
  for (feature in rownames(results_df)) {
    # get contingency categories
    pval <- results_df[feature, pval_colname]
    label <- results_df[feature, label_colname]
    category_vector <- c(category_vector, get_contingency_category(pval, label, threshold))
  }
  df2return <- cbind(results_df, category_vector)
  colnames(df2return)[which(colnames(df2return)=="category_vector")] <- sprintf("category.%s", pval_colname)
  return(df2return)
}

add_DD_performance_to_results <- function(results, annot, FC, threshold, disp_colname, pval_colname) {
  features_in_results <- rownames(results)
  ### add labels to results
  log_FC <- log(FC, 2)
  labels.DD <- ifelse(abs(log(annot[features_in_results, sprintf("%s.S2", disp_colname)]/annot[features_in_results, sprintf("%s.S1", disp_colname)], 2)) > log_FC, 1, 0)
  results <- cbind(results, labels.DD)
  ### add true DD category, i.e. either DD+, DD- or nonDD
  category.DD <- c()
  for (feature in rownames(results)) {
    trueFC.dispersion <- log(annot[feature, sprintf("%s.S2", disp_colname)] / annot[feature, sprintf("%s.S1", disp_colname)], 2)
    category.DD <- c(category.DD, DD_category(trueFC.dispersion, log_FC, "DD"))
  }
  results <- cbind(results, category.DD)
  ### add confusion matrix categories (TP, FP, TN, FN) to results
  results <- add_contingency_category(results, pval_colname, "labels.DD", threshold)
  return(results)
}

performance_stats <- function(results, category_colname) {
  # true positive
  TP <- length(rownames(results[which(results[[category_colname]]=="TP"),]))
  # false negative
  FN <- length(rownames(results[which(results[[category_colname]]=="FN"),]))
  # true negative
  TN <- length(rownames(results[which(results[[category_colname]]=="TN"),]))
  # false positive
  FP <- length(rownames(results[which(results[[category_colname]]=="FP"),]))
  # TPR: True Positive Rate
  TPR <- TP / (TP + FN)
  # TNR: True Negative Rate
  TNR <- TN / (TN + FP)
  # precision (positive predictive value)
  PPV <- TP / (TP + FP)
  # negative predictive value
  NPV <- TN / (TN + FN)
  # false negative rate
  FNR <- FN / (FN + TP)
  if (isFALSE(all.equal(FNR, 1-TPR))) { # need to use all.equal function to compare two floats
    print("Error: FNR is not equal to 1-TPR.")
  }
  # false positive rate
  FPR <- FP / (FP + TN)
  if (isFALSE(all.equal(FPR, 1-TNR))) {
    print("Error: FPR is not equal to 1-TNR.")
  }
  # FDR: False Discovery Rate
  FDR <- FP / (FP + TP)
  if (isFALSE(all.equal(FDR, 1-PPV))) {
    print("Error: FDR is not equal to 1-PPV.")
  }
  # FOR: False Omission Rate
  FOR <- FN / (FN + TN)
  if (isFALSE(all.equal(FOR, 1-NPV))) {
    print("Error: FOR is not equal to 1-NPV.")
  }
  # ACC: accuracy
  ACC <- (TP + TN) / ((TP + TN) + (FP + FN))
  df2return <- data.frame(stat=c("TP", "FN", "TN", "FP", "TPR", "TNR", "PPV", "NPV", "FNR", "FPR", "FDR", "FOR", "ACC"), value=c(TP, FN, TN, FP, TPR, TNR, PPV, NPV, FNR, FPR, FDR, FOR, ACC))
  return(df2return)
}

auc_stats <- function(results, pval_col, label_col, annot, mean_colname, DE_FC) {
  all_auc <- c()
  DE_cat_vector <- c()
  
  # all features
  pred <- prediction(1-results[[pval_col]], results[[label_col]], label.ordering=c(0,1))
  auc <- performance(pred, measure="auc")
  auc.val <- auc@y.values[[1]]
  all_auc <- c(all_auc, auc.val)
  DE_cat_vector <- c(DE_cat_vector, "all")
  
  # per DE categories
  ## add categories according to DE
  features_in_results <- rownames(results)
  log_FC <- log(DE_FC, 2)
  labels.DE <- ifelse(abs(log(annot[features_in_results, sprintf("%s.S2", mean_colname)]/annot[features_in_results, sprintf("%s.S1", mean_colname)], 2)) > log_FC, "DE", "nonDE")
  results <- cbind(results, labels.DE)
  for (one_cat in levels(results[["labels.DE"]])) {
    sub_df <- results[which(results[["labels.DE"]] == one_cat),]
    sub_df <- droplevels(sub_df)
    pred <- prediction(1-sub_df[[pval_col]], sub_df[[label_col]], label.ordering=c(0,1))
    auc <- performance(pred, measure="auc")
    auc.val <- auc@y.values[[1]]
    all_auc <- c(all_auc, auc.val)
    DE_cat_vector <- c(DE_cat_vector, one_cat)
  }
  
  auc_df <- data.frame(category=DE_cat_vector, AUC=all_auc)
  return(auc_df)
}
