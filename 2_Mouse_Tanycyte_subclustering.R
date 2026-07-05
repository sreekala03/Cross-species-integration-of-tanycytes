##### Mouse Tanycyte Subclustering from MEPV object #####

setwd("/Single_cell_manuscript")

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
  library(scProportionTest)
})

#-----------------------------
# Paths
#-----------------------------
main_rds <- "MEPV/MEPV_cells.rds"
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
Idents(Tanycytes_sub) <-  "seurat_clusters"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(5)))              <- "mT.1"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(3)))              <- "mT.2"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(6)))              <- "mT.3"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(1)))              <- "mT.4"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(0)))              <- "mT.5"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(4)))              <- "mT.6"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(2)))              <- "mT.7"
Idents(object = Tanycytes_sub, cells = WhichCells(Tanycytes_sub, ident = c(7)))              <- "mT.8"

Tanycytes_sub$Tany.cell.state <- Idents(object = Tanycytes_sub)
plot1 = DimPlot(Tanycytes_sub, reduction = "humap") & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = F) 

#Tanycytes_sub <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part1_Mouse_tany/Tanycytes_sub.rds")

Idents(Tanycytes_sub) <- factor(
  Idents(Tanycytes_sub),
  levels = c(
    "mT.1","mT.2","mT.3","mT.4",
    "mT.5","mT.6","mT.7","mT.8"
  )
)

Idents(Tanycytes_sub) <- "Tany.cell.state"

# Check cell counts per mT state
print(table(Tanycytes_sub$Tany.cell.state))

# Labelled UMAP
plot_mt <- DimPlot(Tanycytes_sub, reduction = "humap", group.by = "Tany.cell.state") & NoLegend()
plot_mt_labelled <- LabelClusters(plot_mt, id = "ident", size = 4, repel = FALSE)

pdf(file.path(out_dir, "Dimplot_Tanycyte_mT_labelled.pdf"), width = 5.7, height = 4.9)
print(plot_mt_labelled)
dev.off()

# dittoSeq UMAP
pdf(file.path(out_dir, "Fig1d_DittoDimplot_Tanycyte_mT.pdf"), width = 5.7, height = 4.9)
print(dittoDimPlot(Tanycytes_sub, var = "Tany.cell.state", reduction.use = "humap", opacity = 0.7))
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

mouse_genes_to_plot2 <- c(
  "Gpr37l1","Ntsr2","Rfx4","Sema6d","Slc1a2",
  "Shisa9","Rgs20","Rspo3", "Tgfb2", "Flt1","Nr2f2","Adgb","Vwa5b1",
  "Dnah11","Tppp3","Ctxn1","Hipk1",
  "Vcan","Crym","Ephb1","Rorb",
  "Fgf10","Tox2","Mob3b","Frzb","Adamtsl1","Gria3","Dio2","Gpc3","Deptor",
  "Col25a1","Mest",
  "A2m","B2m","Cd44","Ifi44","Ifitm2","Ifitm3","Ifi27", "Gpx1","Cited1",
  "H2-K1",      
  "Cd74",
  "Irf7",
  "Cxcl10"       
  
  dot_mouse_markers = DotPlot(object = Tanycytes_sub, 
        features = mouse_genes_to_plot2,
        assay = "SCT",
        scale = T,
        cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression",
                                 order = 1)) +  RotatedAxis() 

dot_mouse
ggsave(
  filename = file.path(out_dir, "Dotplot_mouse_markers.pdf"),
  plot = dot_mouse_markers,
  width = 13.7,
  height = 3.2,
  units = "in"
)

#-----------------------------
# Highlight mT.8
#-----------------------------
pdf(file.path(out_dir, "Fig1h_Highlight_mT8.pdf"), width = 6, height = 4.8)
print(
  DimPlot(
    Tanycytes_sub,
    reduction = "humap",
    cells.highlight = WhichCells(Tanycytes_sub, idents = "mT.8"),
    sizes.highlight = 1.5
  )
)
dev.off()

#Feature plot
pdf(file.path(out_dir, "Fig.1i_Ifitm3.pdf"), width = 4, height = 4)
print(FeaturePlot(Tanycytes_sub, features = "Ifitm3", reduction = "humap", order = TRUE))
dev.off()

#-----------------------------
# Save final object
#-----------------------------
saveRDS(Tanycytes_sub, "MEPV/Tanycytes_sub.rds")

#-----------------------------
# scProportion permutation tests
#-----------------------------
 
prop_test <- sc_utils(Tanycytes_sub)

reference_condition <- "Fem.Chow.Diest"

test_conditions <- c(
  "Fem.Chow.Proest",
  "Fem.Chow.Est",
  "Male.Chow",
  "Male.HFDS",
  "Fem.HFDS"
)

output_dir <- "MEPV/output/permutation_plots"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

prop_test_results <- list()

for (condition in test_conditions) {
  
  comparison_name <- paste(
    reference_condition,
    "vs",
    condition,
    sep = "_"
  )
  
  test_result <- permutation_test(
    prop_test,
    cluster_identity = "Tany.cell.state",
    sample_1 = condition,
    sample_2 = reference_condition,
    sample_identity = "Sample.name_2"
  )
  
  prop_test_results[[comparison_name]] <- test_result
  
  plot_obj <- permutation_plot(test_result)
  
  pdf(
    file = file.path(
      output_dir,
      paste0(comparison_name, ".pdf")
    ),
    width = 5.2,
    height = 2.3
  )
  
  print(plot_obj)
  dev.off()
}
