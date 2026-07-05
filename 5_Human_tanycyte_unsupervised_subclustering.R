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

####Load human_ME_full
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

##seurat cluster 13 and 14 could be immune/oligodendro contamination
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
#Tanycytes_human_clean <- readRDS("/human_ME_full/Part2_Tany_human_subset/Tanycytes_human_clean.rds")

