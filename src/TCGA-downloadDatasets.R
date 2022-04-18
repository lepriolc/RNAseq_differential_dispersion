library(TCGAbiolinks)
library(SummarizedExperiment)

##############
# Parameters #
##############
# list of datasets
datasets <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-THCA")
# data category and type
data_category <- "Transcriptome Profiling"
data_type <- "Gene Expression Quantification"
# type of workflow
workflow_type <- "Counts"
# phenotypes of interest
phenotypes <- c("TP", "NT")
# path to output directory
output_dir <- "./output/TCGA/00-Data"


############
# Analysis #
############
for (one_dataset in datasets) {
  # get expression data
  ## output directory and basename
  dataset_output_dir <- sprintf("%s/%s", output_dir, one_dataset)
  exp_output_dir <- sprintf("%s/harmonized/%s/%s", dataset_output_dir, gsub(" ", "_", data_category), gsub(" ", "_", data_type))
  dataset_output_basename <- sprintf("%s_%s", one_dataset, workflow_type)
  ## query
  query <- GDCquery(project = one_dataset,
                    data.category = data_category,
                    data.type = data_type,
                    workflow.type = paste("HTSeq - ", workflow_type, sep=""))
  ## download data
  GDCdownload(query, method="client", directory=output_dir)
  ## write query results
  write.csv(query$results[[1]], file=sprintf("%s/%s_query_results.csv", exp_output_dir, dataset_output_basename), quote=FALSE, row.names=FALSE)
  ## get expression matrix
  exp_data <- GDCprepare(query, directory=output_dir)
  exp_matrix <- assay(exp_data, paste("HTSeq - ", workflow_type, sep=""))
  ## get phenotypes
  query_info <- query$results[[1]]
  pheno_df <- data.frame()
  for (pheno in phenotypes) {
    pheno_samples <- TCGAquery_SampleTypes(query_info$cases, pheno)
    one_pheno_df <- data.frame(barcode=pheno_samples, phenotype=rep(pheno, length(pheno_samples)))
    pheno_df <- rbind(pheno_df, one_pheno_df)
  }
  
  # get batch information
  data_category <- "Biospecimen"
  biospecimen_info <- "aliquot"
  ## output directory
  biospecimen_output_dir <- sprintf("%s/harmonized/Biospecimen/Biospecimen_Supplement", dataset_output_dir)
  if (! dir.exists(biospecimen_output_dir)) {
    dir.create(biospecimen_output_dir, recursive=TRUE, mode="0775")
  }
  ## query
  query_biospecimen <- GDCquery(project=one_dataset, data.category=data_category, data.type = "Biospecimen Supplement", file.type="xml")
  ## download data
  GDCdownload(query_biospecimen, method="client", directory=output_dir)
  prepare_clinic_info <- GDCprepare_clinic(query_biospecimen, clinical.info=biospecimen_info, directory=output_dir)
  ## prepare_clinic_info may contain duplicated lines, remove duplicates
  prepare_clinic_info <- prepare_clinic_info[! duplicated(prepare_clinic_info),]
  write.csv(prepare_clinic_info, file=sprintf("%s/%s_%s_%s.csv", dataset_output_dir, one_dataset, data_category, biospecimen_info))
  
  # create metadata table
  ## get metadata in Biospecimen aliquot
  mRNA_barcodes <- as.character(pheno_df$barcode)
  exp_aliquot_df <- prepare_clinic_info[prepare_clinic_info$bcr_aliquot_barcode %in% mRNA_barcodes,]
  metadata_df <- exp_aliquot_df[,c("bcr_aliquot_barcode", "bcr_patient_barcode", "center_id", "plate_id", "plate_row", "plate_column", "source_center")]
  metadata_df$plate_column <- as.factor(metadata_df$plate_column)
  metadata_df$plate_id <- as.factor(metadata_df$plate_id)
  metadata_df$source_center <- as.factor(metadata_df$source_center)
  metadata_df <- droplevels(metadata_df)
  ## get TP and NT phenotypes
  metadata_pheno_df <- merge(pheno_df, metadata_df, by.x="barcode", by.y="bcr_aliquot_barcode")
  ## order by phenotype
  metadata_pheno_df <- metadata_pheno_df[order(metadata_pheno_df$phenotype),]
  rownames(metadata_pheno_df) <- metadata_pheno_df$barcode
  ## rename batch column
  colnames(metadata_pheno_df) <- sub("plate_id", "batch", colnames(metadata_pheno_df))
  
  # filter expression matrix and metadata table: only keep samples from patients with both TP and NT samples
  patient_vector <- as.character(metadata_pheno_df$bcr_patient_barcode)
  ## patients with both TP and NT samples
  patients_more_1_sample <- names(which(table(patient_vector) > 1))
  patient2keep <- unlist(lapply(patients_more_1_sample, function(x) {
    patient_samples_pheno <- metadata_pheno_df[which(metadata_pheno_df$bcr_patient_barcode==x), "phenotype"]
    if ("TP" %in% patient_samples_pheno & "NT" %in% patient_samples_pheno) {
      val2return <- TRUE
    } else {
      val2return <- FALSE
    }
    return(val2return)
  }))
  patients_TP_NT_sample <- patients_more_1_sample[patient2keep]
  ## reduce metadata table and expression matrix
  metadata_pheno_df <- metadata_pheno_df[which(metadata_pheno_df$bcr_patient_barcode %in% patients_TP_NT_sample),]
  metadata_pheno_df <- droplevels(metadata_pheno_df)
  keep_samples <- as.character(metadata_pheno_df$barcode)
  exp_matrix <- exp_matrix[, keep_samples]
  ## reorder expression matrix columns and metadata table lines
  pheno_1 <- levels(metadata_pheno_df$phenotype)[1]
  pheno_2 <- levels(metadata_pheno_df$phenotype)[2]
  pheno_1_samples <- as.character(metadata_pheno_df[which(metadata_pheno_df$phenotype==pheno_1),"barcode"])
  pheno_2_samples <- as.character(metadata_pheno_df[which(metadata_pheno_df$phenotype==pheno_2),"barcode"])
  exp_matrix <- exp_matrix[, c(pheno_1_samples, pheno_2_samples)]
  metadata_pheno_df <- metadata_pheno_df[c(pheno_1_samples, pheno_2_samples),]
  metadata_pheno_df <- droplevels(metadata_pheno_df)
  
  # save expression matrix and metadata table in a RData file
  save(exp_matrix, metadata_pheno_df, file=sprintf("%s/%s_exp_pheno-patients_TP_NT_samples.RData", dataset_output_dir, dataset_output_basename))
}

