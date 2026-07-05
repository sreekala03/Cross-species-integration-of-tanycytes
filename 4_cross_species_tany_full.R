options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)

library(Seurat)
library(Matrix)
library(dplyr)
library(SeuratDisk)

############################
## 1. Load MGI/HMD ortholog table
############################

hmd <- read.delim(
  "HMD_HumanPhenotype.rpt.txt",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(hmd) <- c(
  "Human_Marker_Symbol",
  "Human_Entrez_ID",
  "Mouse_Marker_Symbol",
  "MGI_Accession_ID",
  "Mammalian_Phenotype_ID"
)

orthologs <- data.frame(
  mouse = hmd$Mouse_Marker_Symbol,
  human = hmd$Human_Marker_Symbol,
  stringsAsFactors = FALSE
)

orthologs <- orthologs[
  orthologs$mouse != "" &
    orthologs$human != "" &
    !is.na(orthologs$mouse) &
    !is.na(orthologs$human),
]

# remove repeated phenotype rows
orthologs <- unique(orthologs)

# strict 1:1
orthologs <- orthologs[
  !duplicated(orthologs$mouse) &
    !duplicated(orthologs$human),
]

message("Strict 1:1 orthologs: ", nrow(orthologs))
#Strict 1:1 orthologs: 18079
############################
## 2. Conversion function
############################

convert_species_seurat_offline <- function(seurat_object,
                                           orthologs,
                                           assay = "RNA") {
  mat <- GetAssayData(seurat_object, assay = assay, layer = "counts")
  genes <- rownames(mat)
  cells <- colnames(mat)
  
  meta <- seurat_object@meta.data
  meta <- meta[cells, , drop = FALSE]
  rownames(meta) <- cells
  
  stopifnot(identical(colnames(mat), rownames(meta)))
  
  mouse_overlap <- sum(genes %in% orthologs$mouse)
  human_overlap <- sum(genes %in% orthologs$human)
  
  message("Mouse overlap: ", mouse_overlap)
  message("Human overlap: ", human_overlap)
  
  if (mouse_overlap > human_overlap) {
    message("Converting mouse → human")
    
    conv <- orthologs[orthologs$mouse %in% genes, , drop = FALSE]
    
    mat <- mat[conv$mouse, , drop = FALSE]
    rownames(mat) <- conv$human
    
  } else {
    stop("Input does not look like mouse gene space. Refusing to convert.")
  }
  
  mat <- rowsum(as.matrix(mat), group = rownames(mat))
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  stopifnot(identical(colnames(mat), rownames(meta)))
  
  obj <- CreateSeuratObject(
    counts = mat,
    meta.data = meta,
    project = "mouse_converted_to_human"
  )
  
  stopifnot(identical(colnames(obj), rownames(obj@meta.data)))
  
  return(obj)
}

############################
## 3. Load mouse and human data
############################

mouse <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part1_Mouse_tany/Tanycytes_sub.rds")
human <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part2_Tany_human_subset/Tanycytes_human_clean.rds")

DefaultAssay(mouse) <- "RNA"
DefaultAssay(human) <- "RNA"

stopifnot(identical(colnames(mouse), rownames(mouse@meta.data)))
stopifnot(identical(colnames(human), rownames(human@meta.data)))

mouse$species <- "mouse"
human$species <- "human"

############################
## 4. Convert mouse genes to human orthologs
############################

mouse_conv <- convert_species_seurat_offline(
  seurat_object = mouse,
  orthologs = orthologs,
  assay = "RNA"
)

# Mouse overlap: 16339
# Human overlap: 15
# Converting mouse → human

mouse_conv$species <- "mouse"

#Add tanyycyte labels
colnames(mouse_conv@meta.data)
mouse_conv$tanycyte_label <- as.character(
  mouse_conv$Level2
)
############################
## 5. Normalize and compute HVGs separately
############################

human <- NormalizeData(human, verbose = FALSE)
mouse_conv <- NormalizeData(mouse_conv, verbose = FALSE)

human <- FindVariableFeatures(
  human,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

mouse_conv <- FindVariableFeatures(
  mouse_conv,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

human_hvg <- VariableFeatures(human)
mouse_hvg <- VariableFeatures(mouse_conv)

shared_hvgs <- intersect(human_hvg, mouse_hvg)

message("Shared HVGs before marker rescue: ", length(shared_hvgs))
#Shared HVGs before marker rescue: 1218

scUtils::writeList_to_JSON(
  shared_hvgs,
  file.path(out_dir, "shared_hvgs.json"))


shared_hvgs <- sort(shared_hvgs)
human_sub <- subset(human, features = shared_hvgs)
mouse_sub <- subset(mouse_conv, features = shared_hvgs)
stopifnot(identical(rownames(human_sub), rownames(mouse_sub)))

##Omit computing shared hvgs amd merge objects with full features for post scVI analysis 
############################
## 9. Prefix metadata columns
############################

prefix_metadata <- function(seu, prefix, keep = c("species")) {
  md <- seu@meta.data
  rename_cols <- setdiff(colnames(md), keep)
  colnames(md)[colnames(md) %in% rename_cols] <- paste0(prefix, "_", rename_cols)
  seu@meta.data <- md
  return(seu)
}

human_sub <- prefix_metadata(human_sub, "human")
mouse_sub <- prefix_metadata(mouse_sub, "mouse")

############################
## 10. Rename cells
############################

human_sub <- RenameCells(human_sub, add.cell.id = "human")
mouse_sub <- RenameCells(mouse_sub, add.cell.id = "mouse")

rownames(human_sub@meta.data) <- colnames(human_sub)
rownames(mouse_sub@meta.data) <- colnames(mouse_sub)

stopifnot(identical(colnames(human_sub), rownames(human_sub@meta.data)))
stopifnot(identical(colnames(mouse_sub), rownames(mouse_sub@meta.data)))

############################
## 11. Merge
############################

merged <- merge(
  x = human_sub,
  y = mouse_sub,
  merge.data = FALSE
)

stopifnot(identical(colnames(merged), rownames(merged@meta.data)))

# Join layers first
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

# Convert Assay5 -> old-style Assay
merged[["RNA"]] <- as(merged[["RNA"]], Class = "Assay")

DefaultAssay(merged) <- "RNA"

#Coalesce sample ids mouse and human
merged$sample_id <- coalesce(
  merged$human_sample,
  merged$mouse_Sample.name
)
table(merged@meta.data[["sample_id"]])

#Coalesce cluster ids mouse and human
merged$cluster_id <- coalesce(
  merged$human_seurat_clusters,
  merged$mouse_Tany.rename
)
table(merged@meta.data[["cluster_id"]])

#Add tanycyte labels for mouse with human unknown
merged$mouse_tanycyte_label[
  merged$mouse_tanycyte_label == "ARH tanycytes"
] <- "B1_tanycytes"

merged$mouse_tanycyte_label[
  merged$mouse_tanycyte_label == "DMH tanycytes"
] <- "A1_tanycytes"

merged$mouse_tanycyte_label[
  merged$mouse_tanycyte_label == "ME tanycytes"
] <- "B2_tanycytes"

merged$mouse_tanycyte_label[
  merged$mouse_tanycyte_label == "VMH/ARH tanycytes"
] <- "A2_tanycytes"
table(merged$mouse_tanycyte_label)

#Transfer label Leve1 1 classification
merged$transfer_label <- "Unknown"
mouse_cells <- merged$species == "mouse"
merged$transfer_label[mouse_cells] <-
  merged$mouse_tanycyte_label[mouse_cells]
table(
  merged$species,
  merged$transfer_label
)

#Label transfer with mT classification
merged$transfer_label_2 <- "Unknown"
mouse_cells <- merged$species == "mouse"
merged$transfer_label_2[mouse_cells] <-
  merged$mouse_Tany.rename[mouse_cells]
table(
  merged$species,
  merged$transfer_label_2
)

############################
## 13. Save Seurat object
############################

library(SeuratDisk)
saveRDS(
  merged,
  file = "human_ME_full/Part3_cross_species/before_scVI/tanycytes_cross_species_merged_full.rds")

write.csv(
  merged@meta.data,
  file = "human_ME_full/Part3_cross_species/before_scVI/tanycytes_cross_species_rebuilt_metadata.csv")


############################
## 14. Export h5ad for scVI
############################

# SeuratDisk works better with old-style Assay
merged[["RNA"]] <- as(merged[["RNA"]], Class = "Assay")

SeuratDisk::SaveH5Seurat(
  merged,
  filename = "human_ME_full/tanycytes_cross_species_full.h5seurat",
  overwrite = TRUE
)

SeuratDisk::Convert(
  "human_ME_full/tanycytes_cross_species_full.h5seurat",
  dest = "h5ad",
  assay = "RNA",
  overwrite = TRUE
)
