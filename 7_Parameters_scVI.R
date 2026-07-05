#!/usr/bin/env Rscript

#-----------------------------
# Make scVI/scANVI parameter JSON
#-----------------------------

library(jsonlite)

#-----------------------------
# Paths
#-----------------------------

base_path <- path.expand("~/snRNA_seq_integ/cross_species_tanycytes_LT2")

param_file <- file.path(
  base_path,
  "parameters_cross_species_tanycytes_scvi_LT2.json"
)

merged_file <- file.path(
  base_path,
  "tanycytes_cross_species_LT2.rds"
)

h5ad_file <- file.path(
  base_path,
  "tanycytes_cross_species_LT2.h5ad"
)

feature_set_file <- file.path(
  base_path,
  "shared_hvgs.json"
)

#-----------------------------
# Load feature set
#-----------------------------

shared_hvgs <- unlist(read_json(feature_set_file))

#-----------------------------
# Build parameter list
#-----------------------------

param_list <- list(
  # Paths
  harmonization_folder_path = base_path,
  merged_file = merged_file,
  data_filepath_full = h5ad_file,
  feature_set_file = feature_set_file,
  new_name_suffix = "cross_species_tanycytes_scvi_species_batch",
  
  # General settings
  job_id = "tanycytes_cross_species_LT2",
  n_cores = 8,
  global_seed = 123456,
  
  # Metadata columns
  id_column = "cell_id",
  sample_column = "sample_id",
  batch_var = "species",
  label_key = "transfer_label_2",
  unknown_label = "Unknown",
  
  # Data settings
  feature_set_size = length(shared_hvgs),
  assay_name = "RNA",
  integration_name = "scvi",
  
  # scVI settings
  categorical_covariates = character(0),
  continuous_covariates = character(0),
  n_layers = 2,
  n_latent = 80,
  n_hidden = 256,
  dropout_rate = 0.1,
  max_epochs = 500,
  early_stopping = FALSE,
  dispersion = "gene",
  gene_likelihood = "zinb",
  use_cuda = FALSE,
  
  # Graph / UMAP settings
  k_param = 30,
  dist_type = "cosine",
  
  # scANVI label transfer settings
  scanvi_max_epochs = 500,
  preserve_reference_species = "mouse",
  query_species = "human",
  final_label_column = "scanvi_final_label",
  prediction_column = "scanvi_model_prediction",
  confidence_column = "scanvi_confidence"
)

#-----------------------------
# Sanity checks
#-----------------------------

required_files <- c(
  merged_file,
  h5ad_file,
  feature_set_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Feature set size: ", length(shared_hvgs))
message("Writing parameter file to: ", param_file)

#-----------------------------
# Save JSON
#-----------------------------

write_json(
  param_list,
  path = param_file,
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Done.")
