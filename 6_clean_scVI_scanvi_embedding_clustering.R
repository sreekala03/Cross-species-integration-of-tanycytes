library(Seurat)
library(ggplot2)

set.seed(1234)

out_dir <- "human_ME_full/Part3_cross_species/after_scvi"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------
# Load Seurat object and scANVI outputs
#-----------------------------

merged <- readRDS(
  "human_ME_full/Part3_cross_species/before_scvi/merged_full_features_for_postscVI.rds"
)

scanvi_meta <- read.csv(
  file.path(out_dir, "scanvi_metadata_full.csv"),
  row.names = 1
)

scanvi_latent <- read.delim(
  file.path(out_dir, "scanvi_latent_full.txt"),
  row.names = 1,
  check.names = FALSE
)

#-----------------------------
# Match cells
#-----------------------------

common_cells <- intersect(colnames(merged), rownames(scanvi_latent))

merged_scanvi <- subset(
  merged,
  cells = common_cells
)

scanvi_latent <- scanvi_latent[colnames(merged_scanvi), , drop = FALSE]
scanvi_meta <- scanvi_meta[colnames(merged_scanvi), , drop = FALSE]

merged_scanvi$scanvi_final_label <- scanvi_meta$scanvi_final_label
merged_scanvi$scanvi_confidence <- scanvi_meta$scanvi_confidence

#-----------------------------
# Add scANVI latent space to Seurat
#-----------------------------

scanvi_reduction <- CreateDimReducObject(
  embeddings = as.matrix(scanvi_latent),
  assay = "RNA",
  key = "scANVI_"
)

merged_scanvi[["scanvi"]] <- scanvi_reduction

#-----------------------------
# Run UMAP using scANVI latent space
#-----------------------------

merged_scanvi <- RunUMAP(
  merged_scanvi,
  reduction = "scanvi",
  dims = 1:ncol(Embeddings(merged_scanvi, "scanvi")),
  reduction.name = "umap_scanvi",
  reduction.key = "umapSCANVI_",
  n.neighbors = 30,
  metric = "cosine",
  seed.use = 1234
)

#-----------------------------
# Format labels
#-----------------------------

merged_scanvi$scanvi_final_label <- factor(
  merged_scanvi$scanvi_final_label,
  levels = paste0("mT.", 1:8)
)

merged_scanvi$species <- factor(
  merged_scanvi$species,
  levels = c("mouse", "human")
)

merged_scanvi$T_cluster <- gsub(
  "^mT\\.(\\d+)$",
  "T\\1",
  as.character(merged_scanvi$scanvi_final_label)
)

merged_scanvi$T_cluster <- factor(
  merged_scanvi$T_cluster,
  levels = paste0("T", 1:8)
)

#-----------------------------
# Colors
#-----------------------------

T_cols <- c(
  "T1" = "#E69F00",
  "T2" = "#56B4E9",
  "T3" = "#009E73",
  "T4" = "#F0E442",
  "T5" = "#0072B2",
  "T6" = "#D55E00",
  "T7" = "#F781BF",
  "T8" = "#4D4D4D"
)

species_cols <- c(
  "human" = "#C0392B",
  "mouse" = "#2874A6"
)

#-----------------------------
# Plot scANVI labels
#-----------------------------

p1 <- DimPlot(
  merged_scanvi,
  reduction = "umap_scanvi",
  group.by = "T_cluster",
  cols = T_cols,
  pt.size = 0.5,
  alpha = 0.7
)

ggsave(
  filename = file.path(out_dir, "UMAP_merged_scANVI_T_clusters.pdf"),
  plot = p1,
  width = 4.8,
  height = 3.5
)

#-----------------------------
# Plot scANVI labels split by species
#-----------------------------

p2 <- DimPlot(
  merged_scanvi,
  reduction = "umap_scanvi",
  group.by = "T_cluster",
  split.by = "species",
  cols = T_cols,
  pt.size = 0.5,
  alpha = 0.7
)

ggsave(
  filename = file.path(out_dir, "UMAP_merged_scANVI_T_clusters_species_split.pdf"),
  plot = p2,
  width = 9.5,
  height = 4.8
)

#-----------------------------
# Plot species composition
#-----------------------------

p3 <- DimPlot(
  merged_scanvi,
  reduction = "umap_scanvi",
  group.by = "species",
  cols = species_cols,
  pt.size = 0.5,
  alpha = 0.7
)

ggsave(
  filename = file.path(out_dir, "UMAP_merged_scANVI_species.pdf"),
  plot = p3,
  width = 4.8,
  height = 3.5
)

#-----------------------------
# Save object
#-----------------------------

saveRDS(
  merged_scanvi,
  file = file.path(out_dir, "merged_scanvi_with_umap.rds")
)