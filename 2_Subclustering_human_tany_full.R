########Tanycyte subclustering from human_ME_full####

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)
library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(readr)
library(pheatmap)
library(RColorBrewer)

out_dir <- "human_ME_full/Part2_Tany_human_subset"

####Lodad human_ME_full
human_ME_full_obj <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/human_ME_full_obj.rds")

Idents(human_ME_full_obj) <- "Level1"
Tanycytes_human <-subset(human_ME_full_obj, idents=c("Tanycytes"))

# Use RNA assay as input for SCTransform
DefaultAssay(Tanycytes_human) <- "RNA"

# 3. SCTransform normalization
Tanycytes_human <- SCTransform(
  Tanycytes_human,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)

# 5. Run PCA
Tanycytes_human <- RunPCA(Tanycytes_human, verbose = FALSE)

#
Tanycytes_human <- RunHarmony(
Tanycytes_human,
group.by.vars = "sample"
)

# 7. Choose dimensions
ElbowPlot(Tanycytes_human, ndims = 40)

Tanycytes_human  <- FindNeighbors(object = Tanycytes_human , reduction = "harmony", dims = 1:36) 
Tanycytes_human  <- FindClusters(object = Tanycytes_human , resolution = 1) 
Tanycytes_human  <- RunUMAP(object = Tanycytes_human, dims = 1:36, reduction = "harmony", metric="euclidean")
DimPlot(Tanycytes_human, label = T, order = T, pt.size = 0.1) 

Markers_tany_human <-FindAllMarkers(Tanycytes_human, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.5)

FeaturePlot(
  Tanycytes_human,
  features = c("RAX", "VIM", "COL25A1", "CRYM", "FGFR1", "PTPRC", "OPALIN", "PLP1"),
  reduction = "umap",
  order = TRUE
)

##Cluster 13 and 14 could be immune/oligodendro contamination
VlnPlot(
  Tanycytes_human,
  features = c("PTPRC", "OPALIN"),
  group.by = "seurat_clusters"
)

FeaturePlot(
  Tanycytes_human,
  features = c(
    "PTPRC",
    "AIF1",
    "TYROBP",
    "LST1",
    "FCER1G",
    "SPI1"
  )
)

FeaturePlot(
  Tanycytes_human,
  features = c(
    "OPALIN",
    "MBP",
    "MOG",
    "PLP1",
    "MOBP",
    "MAG"
  ))
VlnPlot(
  Tanycytes_human,
  features = "nCount_RNA"
)

#########Remove clusters 13 and 14; RPL, RPSL and sex genes########

Tanycytes_human_clean <-subset(Tanycytes_human, idents=c("13", "14"), invert = T)

sex.genes <- c(
  "XIST", "TSIX",
  "RPS4Y1", "KDM5D", "UTY",
  "DDX3Y", "EIF1AY", "ZFY"
)

ribo.genes <- grep(
  "^RP[SL]",
  rownames(Tanycytes_human),
  value = TRUE
)
remove.genes <- unique(
  intersect(
    c(sex.genes, ribo.genes),
    rownames(Tanycytes_human)
  )
)
length(remove.genes)
head(remove.genes)

Tanycytes_human_clean <- subset(
  Tanycytes_human_clean,
  features = setdiff(
    rownames(Tanycytes_human_clean),
    remove.genes
  )
)

DefaultAssay(Tanycytes_human_clean) <- "RNA"
#remove old SCT and reductions
Tanycytes_human_clean[["SCT"]] <- NULL
Tanycytes_human_clean[["pca"]] <- NULL
Tanycytes_human_clean[["harmony"]] <- NULL
Tanycytes_human_clean[["umap"]] <- NULL

Tanycytes_human_clean <- SCTransform(
  Tanycytes_human_clean,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)

DefaultAssay(Tanycytes_human_clean) <- "SCT"
Tanycytes_human_clean <- RunPCA(Tanycytes_human_clean, verbose = FALSE)
Tanycytes_human_clean <- RunHarmony(Tanycytes_human_clean, group.by.vars = "sample")


####Cluster again###
# Choose dimensions
ElbowPlot(Tanycytes_human_clean, ndims = 40)

Tanycytes_human_clean  <- FindNeighbors(object = Tanycytes_human_clean, reduction = "harmony", dims = 1:45) 
Tanycytes_human_clean  <- FindClusters(object = Tanycytes_human_clean , resolution = 1.2) 
Tanycytes_human_clean  <- RunUMAP(object = Tanycytes_human_clean, dims = 1:45, reduction = "harmony", metric="euclidean")
DimPlot(Tanycytes_human_clean, label = T, order = T, pt.size = 0.1) 

Idents(Tanycytes_human_clean) <- "seurat_clusters"
Markers_tany_human <-FindAllMarkers(Tanycytes_human_clean, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.5)


# Filter significant markers
Markers_tany_human_sort <- subset(
  Markers_tany_human,
  p_val_adj < 0.05 & avg_log2FC > 0
)

write.csv(Markers_tany_human_sort, file.path(out_dir, "Markers_human_tany_clean_sort.csv"))

saveRDS(Tanycytes_human_clean, "Tanycytes_human_clean.rds")
#Tanycytes_human_clean <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part2_Tany_human_subset/Tanycytes_human_clean.rds")

# ####Similarity between human clusters########
# library(ade4)
# library(pheatmap)
# 
# top100_human <- Markers_tany_human_sort %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 100)
# 
# # Cluster x marker-gene binary matrix
# marker_tab <- table(
#   top100_human$cluster,
#   top100_human$gene
# )
# 
# class(marker_tab)
# marker_mat <- unclass(marker_tab)
# 
# class(marker_mat)
# dim(marker_mat)
# # Binary distance; method = 1 in ade4 is Jaccard distance
# D_jaccard <- ade4::dist.binary(marker_mat, method = 1)
# 
# hc_jaccard <- hclust(D_jaccard, method = "complete")
# 
# plot(
#   hc_jaccard,
#   main = "Cluster similarity based on marker-gene overlap",
#   xlab = "",
#   sub = ""
# )
# sim_jaccard <- 1 - as.matrix(D_jaccard)
# 
# pheatmap::pheatmap(
#   sim_jaccard,
#   clustering_distance_rows = D_jaccard,
#   clustering_distance_cols = D_jaccard,
#   main = "Jaccard similarity of marker genes"
# )
# 
# DotPlot(object = Tanycytes_human_clean, 
#         features = c("TNR", "SHISA9", "SEMA6D", "SLC1A2", "RFX4", "AQP4", "RGS20", "SLC7A11", "ADGB", "VWA5B1", "DNAH11", "VCAN", "CRYM", "NELL2", 
#                      "EPHB1", "ADAMTSL1", "GRIA3", "COL25A1", "DIO2", "FAM20A", "FRZB", "GPC3", "FGF10", "CYP1B1", "A2M", "CD44", "SFRP2", "IGFBP5",
#                      "SCN7A", "B2M", "IFITM1", "IFITM2", "IFI27", "CITED1", "HLA-C", "CD74", "CXCL14", "HLA-DRA"),
#         assay = "RNA",
#         scale = T,
#         cols = "RdYlBu") +
#   geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
#   theme(text = element_text(size = 10),
#         axis.text.x = element_text(hjust = 1,
#                                    vjust = 0.5,
#                                    size = 11,
#                                    color = "black"),
#         axis.text.y = element_text(size = 11),
#         legend.text = element_text(size=9))+
#   labs(title = "", x = "", y = "") +
#   guides(colour = guide_colorbar(title = "Scaled average expression",
#                                  order = 1)) + coord_flip()
# 
# #####Compare these clusters with mouse before annotation#####
# ###Load mouse object and recalculate markers
# Tanycytes_sub <- readRDS("D:/PROJECTS/Single_cell_Jan_2021/BATCH1/ANALYSIS_High_fat_ME/singlecell_ME_HFD/Reanalysis_part1/Single_cell_manuscript/MEPV/Tanycytes_sub.rds")
# Tany_mouse_markers <- FindAllMarkers(Tanycytes_sub, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.1)
# top100_mouse <- Tany_mouse_markers %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 100)
# 
# # -----------------------------
# # 1. Prepare human markers
# # -----------------------------
# human_markers <- top50_human %>%
#   as.data.frame() %>%
#   filter(p_val_adj < 0.05, avg_log2FC > 1) %>%
#   mutate(
#     species = "Human",
#     cluster = (cluster),
#     gene = toupper(gene)
#   ) %>%
#   select(species, cluster, gene)
# 
# mouse_markers <- top50_mouse %>%
#   as.data.frame() %>%
#   filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
#   mutate(
#     species = "Mouse",
#     cluster = (cluster),
#     gene = toupper(gene)
#   ) %>%
#   select(species, cluster, gene)
# 
# # -----------------------------
# # 3. Keep only common marker genes
# # -----------------------------
# common_genes <- intersect(human_markers$gene, mouse_markers$gene)
# 
# human_common <- human_markers %>%
#   filter(gene %in% common_genes)
# 
# mouse_common <- mouse_markers %>%
#   filter(gene %in% common_genes)
# 
# markers_common <- bind_rows(mouse_common, human_common)
# 
# # write.csv(
# #   markers_common,
# #   paste0(out_dir, "markers_human_mouse_common.csv",
# #          row.names = FALSE
# #   ))
# 
# # -----------------------------
# # 4. Build binary cluster x gene matrix
# # -----------------------------
# marker_mat <- table(markers_common$cluster, markers_common$gene)
# marker_mat <- as.matrix(marker_mat)
# storage.mode(marker_mat) <- "numeric"
# marker_mat[marker_mat > 0] <- 1
# 
# # -----------------------------
# # 5. Jaccard distance and similarity
# # -----------------------------
# D_jaccard <- proxy::dist(
#   marker_mat,
#   method = "Jaccard"
# )
# 
# hc_jaccard <- hclust(
#   D_jaccard,
#   method = "complete"
# )
# 
# sim_jaccard <- 1 - as.matrix(D_jaccard)
# 
# human_clusters <- grep("^Human_", rownames(sim_jaccard), value = TRUE)
# mouse_clusters <- grep("^Mouse_", colnames(sim_jaccard), value = TRUE)
# 
# hm_sim <- sim_jaccard[human_clusters, mouse_clusters]
# 
# # -----------------------------
# # 6. Dendrogram
# # -----------------------------
# plot(
#   hc_jaccard,
#   main = "Mouse–Human tanycyte cluster similarity",
#   xlab = "",
#   sub = ""
# )
# 
# 
