##### Tanycyte Subclustering #####
## Cleaned script: keeps only mT.1–mT.8 annotation

setwd("D:/PROJECTS/Single_cell_Jan_2021/BATCH1/ANALYSIS_High_fat_ME/singlecell_ME_HFD/Reanalysis_part1/Single_cell_manuscript")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(dittoSeq)
  library(scCustomize)
})

#-----------------------------
# Paths
#-----------------------------
main_rds <- "MEPV/MEPV_cells_v3.rds"
out_dir  <- "MEPV/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------
# Load main object
#-----------------------------
MEPV_cells <- readRDS(main_rds)
print(MEPV_cells)

#-----------------------------
# Subset tanycytes from main clusters
#-----------------------------
Idents(MEPV_cells) <- "seurat_clusters"

Tanycytes <- subset(
  MEPV_cells,
  idents = c(7, 2, 9, 13)
)

#-----------------------------
# First-pass tanycyte clustering
#-----------------------------
Tanycytes <- FindNeighbors(Tanycytes, dims = 1:28, reduction = "harmony")
Tanycytes <- FindClusters(Tanycytes, resolution = 0.6)
Tanycytes <- RunUMAP(
  Tanycytes,
  dims = 1:28,
  reduction = "harmony",
  reduction.name = "humap",
  reduction.key = "hUMAP_",
  metric = "euclidean"
)

pdf(file.path(out_dir, "Fig0_DimPlot_Tany.pdf"), width = 6.2, height = 5.4)
print(DimPlot(Tanycytes, reduction = "humap", label = TRUE))
dev.off()

# First-pass markers
Tanycytes <- PrepSCTFindMarkers(Tanycytes, assay = "SCT", verbose = TRUE)
Tanycyte_markers <- FindAllMarkers(Tanycytes, only.pos = TRUE)
write.csv(Tanycyte_markers, file.path(out_dir, "Tanycyte_cluster_markers.csv"), row.names = FALSE)

saveRDS(Tanycytes, "MEPV/Tanycytes.rds")

#-----------------------------
# Remove cluster 2: Plp1 / cytoplasmic RNA cluster
#-----------------------------
Idents(Tanycytes) <- "seurat_clusters"

Tanycytes_sub <- subset(
  Tanycytes,
  idents = 2,
  invert = TRUE
)

#-----------------------------
# Recluster filtered tanycytes
#-----------------------------
Tanycytes_sub <- FindNeighbors(Tanycytes_sub, dims = 1:28, reduction = "harmony")
Tanycytes_sub <- FindClusters(Tanycytes_sub, resolution = 0.6)
Tanycytes_sub <- RunUMAP(
  Tanycytes_sub,
  dims = 1:28,
  reduction = "harmony",
  reduction.name = "humap",
  reduction.key = "hUMAP_",
  metric = "euclidean"
)

pdf(file.path(out_dir, "Fig0_DimPlot_Tany_sub_clusters.pdf"), width = 6.2, height = 5.4)
print(DimPlot(Tanycytes_sub, reduction = "humap", label = TRUE))
dev.off()

#-----------------------------
# Annotate tanycyte states as mT.1–mT.8
# Original cluster-to-label mapping:
# 5 -> mT.1
# 3 -> mT.2
# 6 -> mT.3
# 1 -> mT.4
# 0 -> mT.5
# 4 -> mT.6
# 2 -> mT.7
# 7 -> mT.8
#-----------------------------
cluster_to_mt <- c(
  "5" = "mT.1",
  "3" = "mT.2",
  "6" = "mT.3",
  "1" = "mT.4",
  "0" = "mT.5",
  "4" = "mT.6",
  "2" = "mT.7",
  "7" = "mT.8"
)

Tanycytes_sub$Tany.rename <- unname(cluster_to_mt[as.character(Tanycytes_sub$seurat_clusters)])
Tanycytes_sub$Tany.rename <- factor(
  Tanycytes_sub$Tany.rename,
  levels = paste0("mT.", 1:8)
)

Idents(Tanycytes_sub) <- "Tany.rename"

# Check cell counts per mT state
print(table(Tanycytes_sub$Tany.rename))

# Labelled UMAP
plot_mt <- DimPlot(Tanycytes_sub, reduction = "humap", group.by = "Tany.rename") & NoLegend()
plot_mt_labelled <- LabelClusters(plot_mt, id = "ident", size = 4, repel = FALSE)

pdf(file.path(out_dir, "Fig4_Dimplot_Tanycyte_mT_labelled.pdf"), width = 5.7, height = 4.9)
print(plot_mt_labelled)
dev.off()

# dittoSeq UMAP
pdf(file.path(out_dir, "Fig4_DittoDimplot_Tanycyte_mT.pdf"), width = 5.7, height = 4.9)
print(dittoDimPlot(Tanycytes_sub, var = "Tany.rename", reduction.use = "humap", opacity = 0.7))
dev.off()

#-----------------------------
# Markers for mT.1–mT.8 states
#-----------------------------
Tanycytes_sub <- PrepSCTFindMarkers(Tanycytes_sub, assay = "SCT", verbose = TRUE)
Tanycyte_sub_markers <- FindAllMarkers(Tanycytes_sub, only.pos = TRUE)
write.csv(
  Tanycyte_sub_markers,
  file.path(out_dir, "Tanycyte_sub_cluster_markers_mT.csv"),
  row.names = FALSE
)

#-----------------------------
# Highlight mT.8
#-----------------------------
pdf(file.path(out_dir, "Fig7_Highlight_mT8.pdf"), width = 6, height = 4.8)
print(
  DimPlot(
    Tanycytes_sub,
    reduction = "humap",
    cells.highlight = WhichCells(Tanycytes_sub, idents = "mT.8"),
    sizes.highlight = 1.5
  )
)
dev.off()

pdf(file.path(out_dir, "Fig7_Highlight_mT8_split_by_sample.pdf"), width = 12, height = 8)
print(
  DimPlot(
    Tanycytes_sub,
    reduction = "humap",
    cells.highlight = WhichCells(Tanycytes_sub, idents = "mT.8"),
    sizes.highlight = 1.5,
    split.by = "Sample.name"
  )
)
dev.off()

#-----------------------------
# Optional feature plot
#-----------------------------
pdf(file.path(out_dir, "CD44.pdf"), width = 4, height = 4)
print(FeaturePlot(Tanycytes_sub, features = "Cd44", reduction = "humap", order = TRUE))
dev.off()

#-----------------------------
# Save final object
#-----------------------------
saveRDS(Tanycytes_sub, "MEPV/Tanycytes_sub_mT.rds")
