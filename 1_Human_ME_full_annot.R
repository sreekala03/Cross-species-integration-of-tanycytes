####Human ME annotation####

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
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#human_ME <- readRDS("D:/PROJECTS/snRNA_humanME/Human_ME/median_eminence_20122024.rds")
human_ME <- readRDS("D:/PROJECTS/snRNseq_integration/human_ME_full/Part2_Tany_human_subset/human_ME_full_obj.rds")
DimPlot(human_ME, label = T, order = T, pt.size = 0.1)
Idents(human_ME) <- "seurat_clusters"
DimPlot(human_ME, label = T, order = T, pt.size = 0.1)

Idents(human_ME) <- "ann_level_3"
DimPlot(human_ME, label = T, order = T, pt.size = 0.1)

################################################LEVEL 1 annotation####################################
Idents(object = human_ME) <- "seurat_clusters"

Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(0)))               <- "Macrophages"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(10, 11)))          <- "Microglia"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(35,16,36)))        <- "Immune.cells"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(28)))              <- "Endothelial"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(38)))              <- "Pericytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(46)))              <- "Myelinating.glia/Schwann"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(48)))              <- "Ant.pituitary"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(5,3)))             <- "Oligodendrocytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(6)))               <- "Oligo.precursors"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(24)))              <- "VLMCs"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(2,4,42,7)))        <- "Astrocytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(23)))              <- "Ependymocytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(21,41,9,8,18)))    <- "Tanycytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(22,34,40,20,31,
                                                                         25,19,39,27,15,26,
                                                                         44,30,12,29,13,14,
                                                                         17,32,1,43,47,37,33,45
)))    <- "Neurons"


human_ME$Level1 <- Idents(object = human_ME)
plot1 = DimPlot(human_ME) & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = F) 

human_ME$Level1 <- factor(human_ME$Level1,levels=c("Neurons","Astrocytes", "Ependymocytes", "Tanycytes", "Microglia", "Macrophages", 
                                                           "Immune.cells", "Endothelial","Pericytes",
                                                           "VLMCs", "Oligo.precursors", "Oligodendrocytes", 
                                                           "Myelinating.glia/Schwann", "Ant.pituitary"))

###Assign colors to clusters####

col.pal <- list()
col.pal$celltype <- c("Astrocytes"="#A6CEE3", "Tanycytes"="#CBD52E", "Ependymocytes"="darkgoldenrod2", 
                      "Neurons"="#CAB2D6", "Oligo.precursors"="#38b000",  "Oligodendrocytes"="#0466c8",
                      "Myelinating.glia/Schwann" = "#B15928","Ant.pituitary" = "#1B9E77", "Endothelial"="red", 
                      "Pericytes"="yellow2",  "VLMCs"="#FB9A99", "Microglia"= "#7570B3", "Macrophages"="#E7298A",
                      "Immune.cells"="#D95F02")



DimPlot(human_ME, group.by = "Level1", cols = col.pal$celltype, label = T) & NoAxes() & NoLegend() 
cell_no_level1 = as.data.frame(table(human_ME$Level1))

pdf("human_ME_full/1_Dim_human_ME_unlabelled.pdf", width = 15, height = 12)
DimPlot(human_ME, group.by = "Level1", cols = (col.pal$celltype), pt.size = 0.5) & NoAxes() & NoLegend()
dev.off()   

Idents(object = human_ME) <- "Level1"
Markers_human_level1<-FindAllMarkers(human_ME, min.pct=0.25, only.pos = TRUE, logfc.threshold = 1)
Markers_human_level1 <- subset(Markers_human_level1, p_val_adj < 0.05)
write.csv(Markers_human_level1, file.path(out_dir, "Table_S2_Markers_human_level1.csv"))

################################################LEVEL 2 annotation####################################

Idents(object = human_ME) <- "seurat_clusters"

Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(0)))         <- "Macrophages"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(10,11)))     <- "Microglia"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(35)))        <- "B.cells"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(16)))        <- "T.cells"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(36)))        <- "Plasma.cells"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(28)))        <- "Endothelial"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(38)))        <- "Pericytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(46)))        <- "Myelinating.glia/Schwann"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(48)))        <- "Ant.pituitary"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(5,3)))       <- "Oligodendrocytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(6)))         <- "Oligodendrocyte.precursors"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(24)))        <- "VLMCs"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(2,4,42,7)))  <- "Astrocytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(23)))        <- "Ependymocytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(21,41))) <- "A1.tanycytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(9))) <- "A2.B1.tanycytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(8))) <- "B2.tanycytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(18))) <- "Immune.Tanycytes"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(22)))    <- "GLUT.RBFOX1"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(40)))  <- "SIM1"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(31)))  <- "LHX1"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(15,39)))  <- "CRHR2"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(1)))  <- "POMC.1"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(14)))  <- "POMC.2"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(33,12)))  <- "AGRP"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(37)))  <- "GAL/GHRH.1"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(19)))  <- "GAL/GHRH.2"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(47)))  <- "TFPI2+/BRS3+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(43)))  <- "AVP+/OXT+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(13)))  <- "TMEM114+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(29)))  <- "PTPRQ+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(17,32)))  <- "TAC3/KISS1"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(34,20,25,30)))  <- "HTR1E+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(27)))  <- "THSD7B+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(44,26)))  <- "CHRNA7+"
Idents(object = human_ME, cells = WhichCells(human_ME, ident = c(45)))  <- "GLUT_SLC17A7"

human_ME$Level2 <- Idents(object = human_ME)
plot2 = DimPlot(human_ME) & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot2, id = "ident", size = 4, repel = F) 
cell_no_level2 = as.data.frame(table(human_ME$Level2))

Idents(object = human_ME) <- "Level2"
Markers_human_level2<-FindAllMarkers(human_ME, min.pct=0.25, only.pos = TRUE, logfc.threshold = 1)
Markers_human_level2 <- subset(Markers_human_level2, p_val_adj < 0.05)
write.csv(Markers_human_level2, "human_ME_full/Markers_human_level2.csv")

celltype_cols <- list()
celltype_cols$cell_type <- c(
  # Immune cells (reds/oranges)
  "Macrophages"               = "#D73027",
  "Microglia"                 = "#FC8D59",
  "B.cells"                   = "#FDAE61",
  "T.cells"                   = "#F46D43",
  "Plasma.cells"              = "#FEE08B",
  
  # Vascular / stromal (greens)
  "Endothelial"               = "#1B9E77",
  "Pericytes"                 = "#66A61E",
  "VLMCs"                     = "#A6D854",
  
  # Glia (blues)
  "Astrocytes"                = "#4575B4",
  "Oligodendrocytes"          = "#313695",
  "Oligodendrocyte.precursors"= "#74ADD1",
  "Myelinating.glia/Schwann"  = "#92C5DE",
  "Ependymocytes"             = "#5E4FA2",
  
  # Tanycytes (teals)
  "A1.tanycytes"     = "#01665E",
  "A2.B1.tanycytes"  = "#35978F",
  "B2.tanycytes"     = "#80CDC1",
  "Immune.Tanycytes" = "#C7EAE5",
  
  # Other non-neuronal
  "Ant.pituitary"             = "#B2ABD2",
  
  # Neurons (purple → pink spectrum)
  "GLUT.RBFOX1"               = "#762A83",
  "SIM1"                      = "#9970AB",
  "LHX1"                      = "#C2A5CF",
  "CRHR2"                     = "#E7D4E8",
  
  # POMC / AGRP axis
  "POMC.1"                    = "#1F78B4",
  "POMC.2"                    = "#6BAED6",
  "AGRP"                      = "#E31A1C",
  
  # Neuroendocrine populations
  "GAL/GHRH.1"                = "#7A0177",
  "GAL/GHRH.2"                = "#C51B8A",
  "TFPI2+/BRS3+"              = "#F768A1",
  "AVP+/OXT+"                 = "#FDE0DD",
  
  # Additional neuronal populations
  "TMEM114+"                  = "#8C510A",
  "PTPRQ+"                    = "#BF812D",
  "TAC3/KISS1"                = "#DFC27D",
  "HTR1E+"                    = "#80B1D3",
  "THSD7B+"                   = "#B3DE69",
  "CHRNA7+"                   = "#FDB462",
  "GLUT_SLC17A7"              = "#BC80BD"
)


human_ME$Level2 <- factor(
  human_ME$Level2,
  levels = c(
    # Astrocytes
    "Astrocytes", 
    
    # Tanycytes
    "A1.tanycytes", "A2.B1.tanycytes", "B2.tanycytes", "Immune.Tanycytes",
    
    # Other glial / non-neuronal populations
    "Ependymocytes",
    "Oligodendrocyte.precursors",
    "Oligodendrocytes",
    "Myelinating.glia/Schwann",
    "Microglia",
    "Macrophages",
    "B.cells",
    "T.cells",
    "Plasma.cells",
    "Endothelial",
    "Pericytes",
    "VLMCs",
    "Ant.pituitary",
    
    # Neuronal populations
    "POMC.1",
    "POMC.2",
    "AGRP",
    "TAC3/KISS1",
    "GAL/GHRH.1",
    "GAL/GHRH.2",
    "CRHR2",
    "SIM1",
    "LHX1",
    "AVP+/OXT+",
    "TFPI2+/BRS3+",
    "TMEM114+",
    "PTPRQ+",
    "HTR1E+",
    "THSD7B+",
    "CHRNA7+",
    "GLUT.RBFOX1",
    "GLUT_SLC17A7"
  )
)

Idents(human_ME) <- human_ME$Level2


DimPlot(human_ME, group.by = "Level2", cols = celltype_cols$cell_type) & NoAxes()  

saveRDS(human_ME, "human_ME_full/human_ME_full_obj.rds")

# pdf("human_ME_full/2_Dim_human_ME_level2_unlabelled.pdf", width = 15, height = 12)
# DimPlot(human_ME, group.by = "Level2", cols = (celltype_cols$cell_type), pt.size = 0.5) & NoAxes() & NoLegend()
# dev.off()
# 
# 
# 
# pdf("human_ME_full/2_Dim_human_ME_level2_labelled.pdf", width = 10.7, height = 5.8)
# DimPlot(human_ME, group.by = "Level2", cols = celltype_cols$cell_type) & NoAxes()  
# dev.off()
# 
# 
# cell.no = as.data.frame(table(human_ME$Level2))
# 
# 
# 




