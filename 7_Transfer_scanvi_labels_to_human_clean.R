#########Generate mouse and human Tanycyte object with scanvi latent embveddings###

library(Seurat)

library(Seurat)

# Load Seurat object
Tanycytes_human_clean <- readRDS(
  "D:/PROJECTS/snRNseq_integration/human_ME_full/Part2_Tany_human_subset/Tanycytes_human_clean.rds"
)

# Load scANVI metadata
scanvi_metadata <- read.csv(
  "D:/PROJECTS/snRNseq_integration/human_ME_full/Part3_cross_species/after_scvi/scanvi_metadata_full.csv",
  row.names = 1,
  check.names = FALSE
)

# Load scANVI latent representation
scanvi_latent <- read.table(
  "D:/PROJECTS/snRNseq_integration/human_ME_full/Part3_cross_species/after_scvi/scanvi_latent_full.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# Create a copy
Tanycytes_human_scanvi <- Tanycytes_human_clean

# Add human_ prefix to all cell names
Tanycytes_human_scanvi <- RenameCells(
  Tanycytes_human_scanvi,
  new.names = paste0("human_", colnames(Tanycytes_human_scanvi))
)

# Check overlap
length(intersect(
  colnames(Tanycytes_human_scanvi),
  rownames(scanvi_metadata)
))

length(intersect(
  colnames(Tanycytes_human_scanvi),
  rownames(scanvi_latent)
))

# Keep only cells present in the Seurat object
scanvi_metadata <- scanvi_metadata[colnames(Tanycytes_human_scanvi), , drop = FALSE]
scanvi_latent <- scanvi_latent[colnames(Tanycytes_human_scanvi), , drop = FALSE]

# Add scANVI metadata
Tanycytes_human_scanvi <- AddMetaData(
  object = Tanycytes_human_scanvi,
  metadata = scanvi_metadata
)

# Add scANVI latent embedding as a dimensional reduction
Tanycytes_human_scanvi[["scanvi"]] <- CreateDimReducObject(
  embeddings = as.matrix(scanvi_latent),
  key = "SCANVI_",
  assay = DefaultAssay(Tanycytes_human_scanvi)
)

# Check result
Tanycytes_human_scanvi
head(Tanycytes_human_scanvi@meta.data)
Embeddings(Tanycytes_human_scanvi, "scanvi")[1:5, 1:5]

# save updated object
saveRDS(
  Tanycytes_human_scanvi,
  "D:/PROJECTS/snRNseq_integration/human_ME_full/Part3_cross_species/after_scvi/Tanycytes_human_with_scanvi.rds"
)

##Add additional cluster resolution
#human_rds <- file.path(out_dir, "Tanycytes_human_with_scanvi.rds")
#human_rds <- human
human <- FindClusters(
  human, resolution = c(0.4,0.6, 0.8)
)

DimPlot(
  human,
  group.by = "SCT_snn_res.0.4",
  label = TRUE
)

DimPlot(
  human,
  group.by = "SCT_snn_res.0.6",
  label = TRUE
)

DimPlot(
  human,
  group.by = "SCT_snn_res.0.8",
  label = TRUE
)

saveRDS(
  human, "D:/PROJECTS/snRNseq_integration/human_ME_full/Part3_cross_species/after_scvi/Tanycytes_human_with_scanvi.rds"
)

Tanycytes_human_with_scanvi <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part3_cross_species/after_scvi/Tanycytes_human_with_scanvi.rds")

Tanycytes_human_with_scanvi <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part3_cross_species/after_scvi/Tanycytes_human_with_scanvi.rds")
label_col <- "scanvi_final_label"
Tanycytes_human_with_scanvi$T_cluster <- gsub(
  "^mT\\.(\\d+)$",
  "T\\1",
  as.character(Tanycytes_human_with_scanvi[[label_col]][,1])
)
Tanycytes_human_with_scanvi$T_cluster <- factor(
  Tanycytes_human_with_scanvi$T_cluster,
  levels = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"))

T_cols <- c(
  "T1"   = "#E69F00",  # gold
  "T2"   = "#56B4E9",   # blue
  "T3"   = "#009E73",  # green
  "T4"   = "#F0E442",  # pink
  "T5"   = "#0072B2",  # sky blue
  "T6"   = "#D55E00",  # brown
  "T7"   = "#F781BF",  # purple
  "T8"   = "#4D4D4D"   # dark grey
)
p_dim =DimPlot(
  Tanycytes_human_with_scanvi,
  group.by = "T_cluster",
  cols = T_cols,
  label = F
)
ggsave(
  filename = file.path(out_dir, "human_harmony_scanvi_umap.pdf"),
  plot = p_dim,
  width = 6.3,
  height = 5.1,
  units = "in"
)
