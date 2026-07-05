####Single cell RNA sequencing mouse MBH clustering#########
##02 11 2023

# -----------------------------
# 0. Setup
# -----------------------------

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(RColorBrewer)
  library(corrplot)
  library(dittoSeq)
})

project_dir <- "/Single_cell_manuscript"
raw_data_dir <- "/data/Raw"
output_dir <- file.path(project_dir, "MEPV", "output")
rds_dir <- file.path(project_dir, "MEPV")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

setwd(project_dir)

# -----------------------------
# 1. Load samples and create Seurat objects
# -----------------------------

sample_ids <- paste0("s", 1:9)

##sex and dissociation bias genes
genes_to_remove <- c(
  "Ehd2", "Espl1", "Jarid1d", "Pnpla4", "Rps4y1", "Xist", "Tsix",
  "Eif2s3y", "Ddx3y", "Uty", "Kdm5d",
  "Fos", "Fosb", "Gstp1", "Egr1", "Jun", "Junb", "Jund",
  "Erh", "Slc25a5", "Pgk1", "Eno1", "Npas4", "Tubb2a",
  "Emc4", "Scg5", "Gm42418"
)

create_sample_object <- function(sample_id, raw_data_dir) {
  message("Processing sample: ", sample_id)

  sample_obj <- CreateSeuratObject(
    counts = Read10X(file.path(raw_data_dir, sample_id)),
    project = "MBH"
  )

  sample_obj$sample <- sample_id

  sample_obj <- PercentageFeatureSet(
    object = sample_obj,
    pattern = "^mt-",
    col.name = "percent.mt"
  )

  sample_obj <- subset(
    sample_obj,
    subset = nCount_RNA < 20000 &
      nFeature_RNA > 500 &
      nFeature_RNA < 4000 &
      percent.mt < 15
  )

  counts <- GetAssayData(sample_obj, assay = "RNA", slot = "counts")
  ribosomal_genes <- rownames(counts) %>% str_subset("^Rp[sl]")

  keep_features <- setdiff(rownames(counts), c(genes_to_remove, ribosomal_genes))
  sample_obj <- subset(sample_obj, features = keep_features)

  return(sample_obj)
}

MEPV_sample_list <- lapply(sample_ids, create_sample_object, raw_data_dir = raw_data_dir)
names(MEPV_sample_list) <- sample_ids

# -----------------------------
# 2. Merge samples
# -----------------------------

MEPV_cells <- Reduce(function(x, y) merge(x, y), MEPV_sample_list)

var_features <- SelectIntegrationFeatures(
  object.list = MEPV_sample_list,
  nfeatures = 3000,
  fvf.nfeatures = 3000
)

MEPV_cells <- SCTransform(MEPV_cells, vars.to.regress = "percent.mt", verbose = FALSE)

MEPV_cells <- RunPCA(MEPV_cells, features = var_features, verbose = FALSE)

pdf(file.path(output_dir, "QC_elbow_plot.pdf"), width = 6, height = 5)
print(ElbowPlot(MEPV_cells))
dev.off()

pdf(file.path(output_dir, "QC_dim_heatmap_47_50.pdf"), width = 10, height = 8)
print(DimHeatmap(MEPV_cells, dims = 47:50, cells = 500, balanced = TRUE))
dev.off()

MEPV_cells <- RunUMAP(MEPV_cells, dims = 1:50, reduction = "pca")

pdf(file.path(output_dir, "QC_merged_umap_by_sample.pdf"), width = 7, height = 6)
print(DimPlot(MEPV_cells, group.by = "sample") + ggtitle("Merged samples"))
dev.off()          
                     
                     
######## 3. Adding metadata and reordering samples #########

MEPV_cells[["Sex"]]<-ifelse(MEPV_cells$sample%in%c("s1","s2","s4"),"Male","Female")
MEPV_cells[["Diet"]]<-sapply(MEPV_cells$sample,function(s)ifelse(s%in%c("s1","s3","s7","s8"),"Chow",
                                                           ifelse(s%in%c("s4","s9"),"HFDR",
                                                                  ifelse(s%in%c("s5"),"HFDS+D","HFDS"))))

MEPV_cells[["Sample.name"]]<-sapply(MEPV_cells$sample,function(s)ifelse(s%in%c("s1"),"Male.Chow",
                                                                  ifelse(s%in%c("s2"),"Male.HFDS",
                                                                         ifelse(s%in%c("s3"),"Fem.Chow.Diest",
                                                                                ifelse(s%in%c("s4"),"Male.HFDR",
                                                                                       ifelse(s%in%c("s5"),"Fem.HFDS+D",
                                                                                              ifelse(s%in%c("s6"),"Fem.HFDS",
                                                                                                     ifelse(s%in%c("s7"),"Fem.Chow.Proest",
                                                                                                            ifelse(s%in%c("s8"),"Fem.Chow.Est",
                                                                                                                   ifelse(s%in%c("s9"),"Fem.HFDR"))))))))))
MEPV_cells[["Sample.name_2"]]<-sapply(MEPV_cells$sample,function(s)ifelse(s%in%c("s1"),"Male.Chow",
                                                                        ifelse(s%in%c("s2"),"Male.HFDS",
                                                                               ifelse(s%in%c("s3"),"Fem.Chow.Diest",
                                                                                      ifelse(s%in%c("s4"),"Male.HFDR",
                                                                                             ifelse(s%in%c("s5", "s6"),"Fem.HFDS",
                                                                                                    ifelse(s%in%c("s7"),"Fem.Chow.Proest",
                                                                                                                  ifelse(s%in%c("s8"),"Fem.Chow.Est",
                                                                                                                         ifelse(s%in%c("s9"),"Fem.HFDR")))))))))




# -----------------------------
# 4. Harmony integration
# -----------------------------

MEPV_cells <- RunHarmony(
  object = MEPV_cells,
  group.by.vars = "Sample.name",
  assay.use = "SCT",
  plot_convergence = TRUE
)

MEPV_cells <- RunUMAP(
  object = MEPV_cells,
  dims = 1:50,
  reduction = "harmony",
  reduction.name = "humap",
  reduction.key = "hUMAP_",
  metric = "euclidean"
)

pdf(file.path(output_dir, "Harmony_umap_by_sample.pdf"), width = 7, height = 6)
print(DimPlot(MEPV_cells, reduction = "humap", group.by = "sample") +
        ggtitle("Harmony integrated, dims 1:50"))
dev.off()

pdf(file.path(output_dir, "Harmony_umap_split_by_sample.pdf"), width = 12, height = 8)
print(DimPlot(MEPV_cells, reduction = "humap", split.by = "sample") +
        ggtitle("Harmony integrated"))
dev.off()

                                      
                                      
# -----------------------------
# 5. Clustering
# -----------------------------

MEPV_cells <- FindNeighbors(
  object = MEPV_cells,
  dims = 1:50,
  reduction = "harmony"
)

MEPV_cells <- FindClusters(
  object = MEPV_cells,
  resolution = 1
)

pdf(file.path(output_dir, "Clusters_harmony_res1_dims1_50.pdf"), width = 8, height = 6)
print(DimPlot(MEPV_cells, reduction = "humap", label = TRUE) +
        ggtitle("Whole integrated dataset, resolution = 1, dims 1:50"))
dev.off()

# -----------------------------
# 6. Cluster markers
# -----------------------------

Idents(MEPV_cells) <- "seurat_clusters"

MEPV_cells <- PrepSCTFindMarkers(MEPV_cells, assay = "SCT", verbose = TRUE)

markers_int <- FindAllMarkers(
  object = MEPV_cells,
  min.pct = 0.25,
  only.pos = TRUE,
  logfc.threshold = 0.4
)

write.csv(
  markers_int,
  file.path(output_dir, "Markers_mouse_MBH_seurat.csv"),
  row.names = FALSE
)

# -----------------------------
# 7. Cell type annotation
# -----------------------------

assign_cluster_labels <- function(seurat_obj, cluster_to_label, output_column) {
  Idents(seurat_obj) <- "seurat_clusters"

  cluster_labels <- as.character(seurat_obj$seurat_clusters)

  for (label in names(cluster_to_label)) {
    cluster_labels[cluster_labels %in% as.character(cluster_to_label[[label]])] <- label
  }

  seurat_obj[[output_column]] <- cluster_labels
  return(seurat_obj)
}

level1_map <- list(
  "VLMC" = c(4, 21, 29),
  "Plvap Endothelial" = 11,
  "Endothelial cells" = 1,
  "Pericytes" = c(14, 24),
  "VSMC" = 23,
  "Ccl5+" = c(22, 31),
  "Microglia" = c(10, 30),
  "CAMs" = c(17, 32),
  "Progenitors" = 15,
  "Differentiating" = 20,
  "Mature" = c(3, 5, 33, 34, 26),
  "Lhb.Npy.Rax+" = 27,
  "Tanycytes" = c(7, 2, 9, 13),
  "Npy+" = 35,
  "Cell membrane projections" = c(8, 16),
  "Avp+" = 19,
  "Oxt+" = 25,
  "Astrocytes" = c(0, 12),
  "Pars.Tuberalis" = c(18, 28),
  "Ependymocytes" = 6
)

level2_map <- list(
  "VLMC.1" = 4,
  "VLMC.2" = 21,
  "Dural fibroblasts" = 29,
  "Plvap Endothelial" = 11,
  "Endothelial cells" = 1,
  "Pericytes" = c(14, 24),
  "VSMC" = 23,
  "Immune cells" = c(22, 31),
  "Microglia" = c(10, 30),
  "CAMs" = c(17, 32),
  "Progenitors" = 15,
  "Differentiating" = 20,
  "Mature" = c(33, 26, 3, 5, 34),
  "Lhb.Npy.Rax+" = 27,
  "DMH tanycytes" = 7,
  "VMH/dmARH tanycytes" = 2,
  "vmARH tanycytes" = 9,
  "ME tanycytes" = 13,
  "Npy+" = 35,
  "Cell membrane projections" = c(8, 16),
  "Avp+" = 19,
  "Oxt+" = 25,
  "Astrocytes" = c(0, 12),
  "Pars.Tuberalis" = 18,
  "Tight junction epithelial" = 28,
  "Ependymocytes" = 6
)

MEPV_cells <- assign_cluster_labels(MEPV_cells, level1_map, "Level1")
MEPV_cells <- assign_cluster_labels(MEPV_cells, level2_map, "Level2")

Idents(MEPV_cells) <- "Level1"

pdf(file.path(output_dir, "Level1_annotation_umap.pdf"), width = 8, height = 6)
plot_level1 <- DimPlot(MEPV_cells, reduction = "humap") &
  theme(legend.text = element_text(size = 6)) &
  NoLegend()
print(LabelClusters(plot_level1, id = "ident", size = 4, repel = TRUE))
dev.off()

Idents(MEPV_cells) <- "Level2"

pdf(file.path(output_dir, "Level2_annotation_umap.pdf"), width = 8, height = 6)
plot_level2 <- DimPlot(MEPV_cells, reduction = "humap") & NoLegend()
print(LabelClusters(plot_level2, id = "ident", size = 4, repel = FALSE))
dev.off()

# Remove Lhb.Npy.Rax+ cells- dissocitaion artefact anterior pituitary 
MEPV_cells <- subset(MEPV_cells, subset = Level2 != "Lhb.Npy.Rax+")

cell_num_v2 <- table(MEPV_cells$Level2, MEPV_cells$Sample.name)
write.csv(cell_num_v2, file.path(output_dir, "Cell_number_v2.csv"))

# -----------------------------
# 9. Reorder identities and define colors
# -----------------------------

level2_order <- c(
  "Astrocytes", "Ependymocytes", "DMH tanycytes", "VMH/dmARH tanycytes",
  "vmARH tanycytes", "ME tanycytes", "Npy+", "Cell membrane projections",
  "Avp+", "Oxt+", "Microglia", "CAMs", "Immune cells",
  "Endothelial cells", "Plvap Endothelial", "Pericytes", "VSMC",
  "VLMC.1", "VLMC.2", "Dural fibroblasts", "Pars.Tuberalis",
  "Tight junction epithelial", "Progenitors", "Differentiating", "Mature"
)

MEPV_cells$Level2 <- factor(MEPV_cells$Level2, levels = level2_order)

celltype_colors <- c(
  "Astrocytes" = "#829399",
  "Ependymocytes" = "darkgoldenrod1",
  "DMH tanycytes" = "#1F78B4",
  "VMH/dmARH tanycytes" = "#B2DF8A",
  "vmARH tanycytes" = "#33A02C",
  "ME tanycytes" = "#FB9A99",
  "Npy+" = "#E31A1C",
  "Oxt+" = "#CAB2D6",
  "Cell membrane projections" = "#6A3D9A",
  "Avp+" = "#CBD52E",
  "Immune cells" = "#B15928",
  "VLMC.1" = "#1B9E77",
  "VLMC.2" = "#D95F02",
  "Dural fibroblasts" = "#CBD52E",
  "Microglia" = "#7570B3",
  "CAMs" = "#E7298A",
  "Pars.Tuberalis" = "#66A61E",
  "Tight junction epithelial" = "#E6AB02",
  "Plvap Endothelial" = "#FF7F00",
  "Endothelial cells" = "#abc4ff",
  "Pericytes" = "yellow2",
  "VSMC" = "tomato1",
  "Mature" = "#0466c8",
  "Progenitors" = "#38b000",
  "Differentiating" = "ivory4"
)

# -----------------------------
# 10. Main UMAP plots
# -----------------------------

pdf(file.path(output_dir, "Fig1_Dimplot_harmony_labelled.pdf"), width = 12.5, height = 8)
print(DimPlot(
  MEPV_cells,
  group.by = "Level2",
  cols = celltype_colors,
  reduction = "humap",
  label = TRUE,
  label.box = TRUE
) & NoLegend() & NoAxes())
dev.off()

saveRDS(MEPV_cells, file.path(rds_dir, "MEPV_cells.rds"))




                                      

