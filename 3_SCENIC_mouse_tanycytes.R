# ============================================================
# SCENIC pipeline for mouse tanycyte subclusters
# ============================================================
# Input:  Tanycytes_sub Seurat object
# Output: SCENIC regulons, AUCell scores, binarized regulon activity,
#         heatmaps and top-regulator tables by sample group and mT state
# ============================================================

# -----------------------------
# 0. Package installation
# -----------------------------
# # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c(
#   "AUCell", "RcisTarget", "GENIE3", "zoo", "mixtools", "rbokeh",
#   "DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne", "doRNG"
# ))
#
# # Optional / platform dependent parallel backend.
# # doMC is not available on Windows.
# BiocManager::install("doMC")
#
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCENIC")
#
# # Required for current RcisTarget feather database reading on some systems.
# options(download.file.method = "libcurl")
# install.packages(
#   "arrow",
#   repos = "https://packagemanager.rstudio.com/all/__linux__/focal/latest"
# )

# -----------------------------
# 1. Setup
# -----------------------------

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1300 * 1024^2)
set.seed(1234)

library(Seurat)
library(dplyr)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(doRNG)
library(foreach)
library(rngtools)
library(doParallel)
library(ComplexHeatmap)
library(reshape2)

project_dir <- "/Tanycyte_scenic"
setwd(project_dir)

dir.create("int", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("output/heatmaps", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Load tanycyte Seurat object
# -----------------------------

Tanycytes_sub <- readRDS(file.path(project_dir, "Tanycytes_sub.rds"))

# Check required metadata columns before downstream summaries.
required_metadata <- c("Sample.name_2", "Tany.cell.state")
missing_metadata <- setdiff(required_metadata, colnames(Tanycytes_sub@meta.data))

if (length(missing_metadata) > 0) {
  warning(
    "The following metadata columns are missing and related summaries may fail: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# -----------------------------
# 3. Initialize SCENIC
# -----------------------------

org <- "mgi"  # mouse; use "hgnc" for human or "dmel" for fly
myDatasetTitle <- "SCENIC_Tanycytes"
dbDir <- file.path(project_dir, "RcisTarget")

# Loads default database names for selected organism.
data(defaultDbNames)
dbs <- defaultDbNames[[org]]

scenicOptions <- initializeScenic(
  org = org,
  dbDir = dbDir,
  dbs = dbs,
  datasetTitle = myDatasetTitle,
  nCores = 10
)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 1234
saveRDS(scenicOptions, file = "int/scenicOptions.rds")

# -----------------------------
# 4. Extract and filter expression matrix
# -----------------------------

exprMat <- as.matrix(GetAssayData(Tanycytes_sub, assay = "RNA", slot = "counts"))
print(dim(exprMat))

# Keep genes expressed in at least 1% of cells and with sufficient total counts.
genesKept <- geneFiltering(
  exprMat,
  scenicOptions = scenicOptions,
  minCountsPerGene = 3 * 0.01 * ncol(exprMat),
  minSamples = 0.01 * ncol(exprMat)
)

exprMat_filtered <- exprMat[genesKept, ]
print(dim(exprMat_filtered))

rm(exprMat)
gc()

# -----------------------------
# 5. Correlation and GENIE3 network inference
# -----------------------------

runCorrelation(exprMat_filtered, scenicOptions)

# SCENIC GENIE3 is typically run on log-transformed expression.
exprMat_filtered <- log2(exprMat_filtered + 1)

# This step can take hours to days depending on dataset size and cores.
runGenie3(exprMat_filtered, scenicOptions)

# -----------------------------
# 6. Build regulons and score cells
# -----------------------------

scenicOptions <- readRDS("int/scenicOptions.rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 1234

# Step 1: infer co-expression modules.
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# Step 2: create regulons using RcisTarget motif enrichment.
scenicOptions <- runSCENIC_2_createRegulons(
  scenicOptions,
  coexMethod = c("top5perTarget", "top10perTarget", "top50perTarget")
)

# Step 3: score regulon activity in individual cells with AUCell.
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

saveRDS(scenicOptions, file = "int/scenicOptions.rds")

# -----------------------------
# 7. Optional: manually adjust AUCell thresholds
# -----------------------------
# This opens a Shiny app. Run only in an interactive R session.

# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
# savedSelections <- shiny::runApp(aucellApp)
#
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds", 1] <- "int/newThresholds.rds"
# saveRDS(newThresholds, file = getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file = "int/scenicOptions.rds")

# -----------------------------
# 8. Binarize regulon activity
# -----------------------------

scenicOptions@settings$devType <- "png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

saveRDS(scenicOptions, file = "int/scenicOptions.rds")

# -----------------------------
# 9. t-SNE on regulon AUC matrix
# -----------------------------

scenicOptions@settings$seed <- 123
nPcs <- c(5, 15, 50)

fileNames_auc <- tsneAUC(
  scenicOptions,
  aucType = "AUC",
  nPcs = nPcs,
  perpl = c(5, 15, 50)
)

fileNames_auc_high_conf <- tsneAUC(
  scenicOptions,
  aucType = "AUC",
  nPcs = nPcs,
  perpl = c(5, 15, 50),
  onlyHighConf = TRUE,
  filePrefix = "int/tSNE_oHC"
)

print(tsneFileName(scenicOptions))

# -----------------------------
# 10. Inspect regulons and motif support
# -----------------------------

regulons <- loadInt(scenicOptions, "regulons")

# Example: inspect specific regulons if present.
regulons[c("Irf7")]

aucell_regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(aucell_regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
write.csv(
  regulonTargetsInfo,
  file = "output/tables/SCENIC_regulonTargetsInfo.csv",
  row.names = FALSE
)

# Example motif inspection for a TF of interest.
# tableSubset <- regulonTargetsInfo[TF == "Stat3" & highConfAnnot == TRUE]
# viewMotifs(tableSubset, options = list(pageLength = 5))
#
# motifEnrichment_selfMotifs_wGenes <- loadInt(
#   scenicOptions,
#   "motifEnrichment_selfMotifs_wGenes"
# )
# tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs == "Irf7"]
# viewMotifs(tableSubset)

# -----------------------------
# 11. Helper function: regulon heatmap and top regulators
# -----------------------------

plot_regulon_activity_by_group <- function(
    seurat_obj,
    regulon_auc,
    group_by,
    output_prefix,
    heatmap_width = 5.9,
    heatmap_height = 15,
    cluster_columns = FALSE
) {
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column not found: ", group_by)
  }

  cell_info <- data.frame(
    cell = colnames(seurat_obj),
    group = seurat_obj@meta.data[[group_by]],
    row.names = colnames(seurat_obj)
  )

  cell_info <- cell_info[!is.na(cell_info$group), , drop = FALSE]

  auc_matrix <- getAUC(regulon_auc)
  common_cells <- intersect(colnames(auc_matrix), rownames(cell_info))

  cell_info <- cell_info[common_cells, , drop = FALSE]
  auc_matrix <- auc_matrix[, common_cells, drop = FALSE]

  regulon_activity <- sapply(
    split(rownames(cell_info), cell_info$group),
    function(cells) rowMeans(auc_matrix[, cells, drop = FALSE])
  )

  regulon_activity_scaled <- t(scale(t(regulon_activity), center = TRUE, scale = TRUE))
  regulon_activity_scaled[is.na(regulon_activity_scaled)] <- 0

  pdf(
    file = file.path("output/heatmaps", paste0(output_prefix, "_heatmap.pdf")),
    width = heatmap_width,
    height = heatmap_height
  )
  ComplexHeatmap::Heatmap(
    regulon_activity_scaled,
    name = "Regulon activity",
    cluster_columns = cluster_columns
  )
  dev.off()

  top_regulators <- reshape2::melt(regulon_activity_scaled)
  colnames(top_regulators) <- c("Regulon", group_by, "RelativeActivity")
  top_regulators <- top_regulators[top_regulators$RelativeActivity > 0, ]

  write.csv(
    top_regulators,
    file = file.path("output/tables", paste0(output_prefix, "_top_regulators.csv")),
    row.names = FALSE
  )

  return(list(
    activity = regulon_activity,
    activity_scaled = regulon_activity_scaled,
    top_regulators = top_regulators
  ))
}

# -----------------------------
# 12. Regulon heatmap by sample condition
# -----------------------------

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

sample_results <- plot_regulon_activity_by_group(
  seurat_obj = Tanycytes_sub,
  regulon_auc = regulonAUC,
  group_by = "Sample.name_2",
  output_prefix = "SCENIC_by_Sample.name",
  cluster_columns = FALSE
)

# -----------------------------
# 13. Regulon heatmap by tanycyte mT state
# -----------------------------
# This uses the mT.1-mT.8 annotation column.

if ("Tany.cell.state" %in% colnames(Tanycytes_sub@meta.data)) {
  mt_results <- plot_regulon_activity_by_group(
    seurat_obj = Tanycytes_sub,
    regulon_auc = regulonAUC,
    group_by = "Tany.cell.state",
    output_prefix = "SCENIC_by_mT_state",
    cluster_columns = TRUE
  )
}

# -----------------------------
# 15. Save workspace
# -----------------------------

save.image("SCENIC_tanycytes.RData")
