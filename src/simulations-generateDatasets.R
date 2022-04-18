library(compcodeR)
library(ggplot2)

##############
# Parameters #
##############
# number of samples per condition
samples_per_condition <- 50
# number of genes
genes <- 10000
# sequencing depth
depth <- 5000000
min_depth_param <- max_depth_param <- 1
# DE genes
DE_fraction <- 0
# DE_fraction <- 0
DE_down_fraction <-  0.5 # fraction of DE with a decrease of mean in the second condition
DE_base_effect <- 1.5
non_DE_runif_FC <- TRUE # generate FC for non DE genes using runif() function
# DD genes
DD_fraction <- 0.5
DD_down_fraction <- 0.5 # fraction of DD with a decrease of dispersion in the second condition
DD_base_effect <- 1.5
DD_effect_param <- 1
non_DD_runif_FC <- FALSE # generate FC for non DD genes using runif() function
# outliers
single_outlier_high_fraction <- 0.1
single_outlier_low_fraction <- 0
random_outlier_high_fraction <- 0
random_outlier_low_fraction <- 0
# number of replicates
replicates <- 10
# path to output directory
output_dir <- "./output/simulations/00-Data"


#############
# Functions #
#############
source("./src/simulations-functions.R")


############
# Analysis #
############
nb_samples_condition_2 <- samples_per_condition
DE_effect_param <- ifelse(DE_fraction==0, 0, get_DE_effect_param(DE_base_effect))
non_DE_runif_FC_name <- ifelse(non_DE_runif_FC, "runif", "1")
non_DD_runif_FC_name <- ifelse(non_DD_runif_FC, "runif", "1")
output_subdir_name <- ifelse(DE_fraction==0, "10-lowly_DE", "00-highly_DE")
output_dir <- sprintf("%s/%s/base_effect_%s/%s", output_dir, output_subdir_name, sub("[.]", "_", DE_base_effect), samples_per_condition)
dataset_basename <- sprintf("%d_%d_%d_%d_DE_%s_%s_%s_%s_DD_%s_%s_%s_%s_SO_%s_%s_RO_%s_%s", genes, samples_per_condition, nb_samples_condition_2, depth, sub("[.]", "_", DE_fraction), sub("[.]", "_", DE_base_effect), sub("[.]", "_", DE_effect_param), non_DE_runif_FC_name, sub("[.]", "_", DD_fraction), sub("[.]", "_", DD_base_effect), sub("[.]", "_", DD_effect_param), non_DD_runif_FC_name, sub("[.]", "_", single_outlier_high_fraction), sub("[.]", "_", single_outlier_low_fraction), sub("[.]", "_", random_outlier_high_fraction), sub("[.]", "_", random_outlier_low_fraction))

for (i in 1:replicates) {
  # replicate output directory
  one_replicate_output_dir <- sprintf("%s/repl_%d", output_dir, i)
  if (! dir.exists(one_replicate_output_dir)) {
    dir.create(one_replicate_output_dir, recursive=TRUE, mode="0775")
  }
  # dataset name
  dataset_name <- sprintf("%s_repl_%d", dataset_basename, i)
  print(dataset_name)
  
  # create a first dataset to have dispersion values from Pickrell and Cheung datasets
  ## same dispersion in the two groups, no outlier counts
  Pickrell_Cheung_dataset_name <- sprintf("Pickrell_Cheung_mean_dispersion_repl_%d", i)
  out <- sprintf("%s/%s", one_replicate_output_dir, Pickrell_Cheung_dataset_name)
  out_rds_file <- sprintf("%s.rds", out)
  out_html_file <- sprintf("%s.html", out)
  effect.size.mean <- create_effect_size_vector(genes, DE_fraction, DE_down_fraction, DE_base_effect, DE_effect_param, non_DE_runif_FC)
  Pickrell_Cheung_dataset <- generateSyntheticData(dataset=Pickrell_Cheung_dataset_name, n.vars=genes, samples.per.cond=samples_per_condition, repl.id=i, 
                                                   seqdepth=depth, minfact=min_depth_param, maxfact=max_depth_param,
                                                   relmeans="auto", effect.size=effect.size.mean, 
                                                   dispersions="auto", between.group.diffdisp=TRUE, 
                                                   filter.threshold.total=0, 
                                                   output.file=out_rds_file)
  summarizeSyntheticDataSet(data.set=out_rds_file, 
                            output.filename=out_html_file)
  write.csv(Pickrell_Cheung_dataset@variable.annotations, file=sprintf("%s/%s_variable_annotations.csv", one_replicate_output_dir, Pickrell_Cheung_dataset_name), quote=FALSE, row.names=TRUE)
  ## get mean and dispersion values in pairs from Pickrell and Cheung datasets for condition 1
  mean_values <- Pickrell_Cheung_dataset@variable.annotations$truemeans.S1[1:genes]
  dispersion_values <- Pickrell_Cheung_dataset@variable.annotations$truedispersions.S1[1:genes]
  
  # simulate dataset
  ## create dispersion values for both conditions
  ### effect size for DE genes
  nb_DE_genes <- floor(genes * DE_fraction)
  nb_DE_down_genes <- floor(nb_DE_genes * DE_down_fraction)
  effect.size.dispersion.DE.down <- create_effect_size_vector(nb_DE_down_genes, DD_fraction, DD_down_fraction, DD_base_effect, DD_effect_param, non_DD_runif_FC)
  nb_DE_up_genes <- nb_DE_genes - nb_DE_down_genes
  effect.size.dispersion.DE.up <- create_effect_size_vector(nb_DE_up_genes, DD_fraction, DD_down_fraction, DD_base_effect, DD_effect_param, non_DD_runif_FC)
  ### effect size for non DE genes
  nb_non_DE_genes <- genes - nb_DE_genes
  nb_non_DE_down_genes <- floor(nb_non_DE_genes/2)
  effect.size.dispersion.nonDE.down <- create_effect_size_vector(nb_non_DE_down_genes, DD_fraction, DD_down_fraction, DD_base_effect, DD_effect_param, non_DD_runif_FC)
  nb_non_DE_up_genes <- nb_non_DE_genes - nb_non_DE_down_genes
  effect.size.dispersion.nonDE.up <- create_effect_size_vector(nb_non_DE_up_genes, DD_fraction, DD_down_fraction, DD_base_effect, DD_effect_param, non_DD_runif_FC)
  ### create dispersion value matrix for all genes
  all.effect.size.dispersion <- c(effect.size.dispersion.DE.down, effect.size.dispersion.DE.up, effect.size.dispersion.nonDE.down, effect.size.dispersion.nonDE.up)
  disp_cond_1 <- dispersion_values
  disp_cond_2 <- all.effect.size.dispersion[1:length(disp_cond_1)] * disp_cond_1
  disp_matrix <- cbind(disp_cond_1, disp_cond_2)
  
  ## generate dataset
  out <- sprintf("%s/%s", one_replicate_output_dir, dataset_name)
  out_rds_file <- sprintf("%s.rds", out)
  out_html_file <- sprintf("%s.html", out)
  sim_dataset <- generateSyntheticData(dataset=dataset_name, n.vars=genes, samples.per.cond=samples_per_condition, repl.id=i,
                                       seqdepth=depth, minfact=min_depth_param, maxfact=max_depth_param, 
                                       relmeans=mean_values, effect.size=effect.size.mean, dispersions=disp_matrix, 
                                       single.outlier.high.prob=single_outlier_high_fraction, single.outlier.low.prob=single_outlier_low_fraction, 
                                       random.outlier.high.prob=random_outlier_high_fraction, random.outlier.low.prob=random_outlier_low_fraction, 
                                       filter.threshold.total=0, filter.threshold.mediancpm=0, fraction.non.overdispersed=0, 
                                       output.file=out_rds_file)
  ### html report
  summarizeSyntheticDataSet(data.set=out_rds_file, 
                            output.filename=out_html_file)
  ### move md files: from working directory to output directory
  dataset_md_file <- sprintf("%s.md", dataset_name)
  file.copy(dataset_md_file, one_replicate_output_dir, copy.date=TRUE)
  unlink(dataset_md_file)
  Pickrell_Cheung_md_file <- sprintf("%s.md", Pickrell_Cheung_dataset_name)
  file.copy(Pickrell_Cheung_md_file, one_replicate_output_dir, copy.date=TRUE)
  unlink(Pickrell_Cheung_md_file)
  ### variable annotations
  variable_annotations <- sim_dataset@variable.annotations
  #### add true log2FC for dispersion
  truelog2foldchanges.dispersion <- log(variable_annotations$truedispersions.S2 / variable_annotations$truedispersions.S1, 2)
  variable_annotations <- cbind(variable_annotations, truelog2foldchanges.dispersion)
  write.csv(variable_annotations, file=sprintf("%s/%s_variable_annotations.csv", one_replicate_output_dir, dataset_name), quote=FALSE, row.names=TRUE)
  #### mean and dispersion scatter plot
  pdf(sprintf("%s/%s_mean_dispersion_FC_scat_plot.pdf", one_replicate_output_dir, dataset_name))
  p <- ggplot(variable_annotations, aes(x=truelog2foldchanges, y=truelog2foldchanges.dispersion)) +
    geom_point() +
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    geom_vline(xintercept=0, linetype="dashed", color="black") +
    labs(title = labs(title = sprintf("True log2FC for mean and dispersion\n#samples condition 1: %d, #samples condition 2: %d\npct DE: %.1f, base effect DE: %.1f, effect parameter DE: %.1f\npct DD: %.1f, base effect DD: %.1f, effect parameter DD: %.1f", samples_per_condition, nb_samples_condition_2, DE_fraction, DE_base_effect, DE_effect_param, DD_fraction, DD_base_effect, DD_effect_param)))
  print(p)
  dev.off()
}

