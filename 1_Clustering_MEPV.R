####Single cell ME PV analysis#########
##02 11 2023

setwd("D:/PROJECTS/Single_cell_Jan_2021/BATCH1/ANALYSIS_High_fat_ME/singlecell_ME_HFD/Reanalysis_part1/Single_cell_manuscript")
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Nebulosa)
library(harmony)
library(readr)
library(venn)
library(dplyr)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(enrichR)
library(rafalib)
library(multtest)
library(metap)
library(RColorBrewer)

MEPV_cells <- readRDS("D:/PROJECTS/Single_cell_Jan_2021/BATCH1/ANALYSIS_High_fat_ME/singlecell_ME_HFD/RDS_files/mito_reg_ribo_stress_v2/samples_mito_reg_ribo_stress_v2.rds")

###### LOAD SAMPLES & CREATE SEURAT OBJECT#######

MEPV_sample_list<-lapply(paste0("s",1:9),function(s){
  sample<-CreateSeuratObject(Read10X(file.path("/data/Raw",s)),project = "MBH")
  sample[["sample"]]<-s
  sample<-PercentageFeatureSet(sample,pattern = "^mt-",col.name = "percent.mt")
  sample<-subset(sample, subset = nCount_RNA < 20000 & nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)
  counts<- GetAssayData(sample, assay="RNA")
  Rpl.genes <- rownames(counts) %>% stringr::str_subset(string = ., pattern = "^Rp[sl]")
  counts<-counts[-(which(rownames(counts) %in% c('Ehd2', 'Espl1', 'Jarid1d', 'Pnpla4',  'Rps4y1', 'Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d', 'Fos', 'Fosb', 'Gstp1', 'Egr1', 'Jun', 'Junb', 'Jund','Erh', 'Slc25a5', 'Pgk1', 'Eno1', 'Npas4', 'Tubb2a', 'Emc4', 'Scg5', Rpl.genes, 'Gm42418'))),]
  sample<-subset(sample, features=rownames(counts))
  return(sample)
})


####### Merge samples ###########

MEPV_cells <-Reduce(merge, MEPV_sample_list)
var_features <- SelectIntegrationFeatures(object.list = MEPV_sample_list, nfeatures = 3000, fvf.nfeatures = 3000)

MEPV_cells <- RunPCA(MEPV_cells,features = var_features)
ElbowPlot(MEPV_cells)

DimHeatmap(MEPV_cells, dims = 47:50, cells = 500, balanced = T)

MEPV_cells <- RunUMAP(object = MEPV_cells, dims = 1:50)
DimPlot(MEPV_cells,group.by="sample")+ggtitle("merged")



######## Adding metadata and reordering samples #########

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


Tanycytes_sub[["Sample.name_2"]]<-sapply(Tanycytes_sub$sample,function(s)ifelse(s%in%c("s1"),"Male.Chow",
                                                                          ifelse(s%in%c("s2"),"Male.HFDS",
                                                                                 ifelse(s%in%c("s3"),"Fem.Chow.Diest",
                                                                                        ifelse(s%in%c("s4"),"Male.HFDR",
                                                                                               ifelse(s%in%c("s5", "s6"),"Fem.HFDS",
                                                                                                      ifelse(s%in%c("s7"),"Fem.Chow.Proest",
                                                                                                             ifelse(s%in%c("s8"),"Fem.Chow.Est",
                                                                                                                    ifelse(s%in%c("s9"),"Fem.HFDR")))))))))





######### INTERGRATION USING HARMONY ################

library(harmony)
MEPV_cells<-RunHarmony(samples,group.by.vars = "Sample.name",assay.use = "SCT", plot_convergence=TRUE)
# Harmony converged after 9 iterations
MEPV_cells <- RunUMAP(object = samples, dims = 1:50, reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric = 'euclidean')
p1 <- DimPlot(samples,reduction="humap",group.by="sample")+ggtitle("Harmony integrated 1:50")
p1

DimPlot(MEPV_cells,reduction="humap",split.by="sample")+ggtitle("Harmony integrated")

########## CLUSTERING ###########

MEPV_cells <- FindNeighbors(object = MEPV_cells, dims = 1:50,reduction = "harmony") 
MEPV_cells <- FindClusters(object = MEPV_cells, resolution = 1) 
MEPV_cells <- RunUMAP(object = MEPV_cells, dims = 1:50,reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric="euclidean")

DimPlot(MEPV_cells, reduction = "humap", label=TRUE)+ggtitle("Whole Integrated dataset, res=1, dim 1:50")

p2 <- DimPlot(MEPV_cells, split.by ="Sample.name", reduction = "humap", label=TRUE)+ggtitle("Whole Integrated dataset, res=0.8")

DimPlot(MEPV_cells, reduction = "humap", label=TRUE)
DimPlot(MEPV_cells_1, reduction = "humap", label=TRUE)

##################### Cluster Markers  #########################

Idents(MEPV_cells) <- "seurat_clusters"

MEPV_cells <- PrepSCTFindMarkers(MEPV_cells, assay = "SCT", verbose = TRUE)
Markers_int<-FindAllMarkers(MEPV_cells, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.4)
write.csv(Markers, "MEPV/output/Markers_int.csv")


############# CELL TYPE ANNOTATION ###############

Idents(object = MEPV_cells) <- "seurat_clusters"

#Level1 annotation 

Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(4, 21, 29)))       <- "VLMC"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(11)))              <- "Plvap Endothelial"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(1)))               <- "Endothelial cells"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(14, 24)))          <- "Pericytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(23)))              <- "VSMC"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(22, 31)))          <- "Ccl5+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(10, 30)))          <- "Microglia" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(17,32)))           <- "CAMs" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(15)))              <- "Progenitors"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(20)))              <- "Differentiating"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(3,5, 33, 34, 26))) <- "Mature"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(27)))              <- "Lhb.Npy.Rax+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(7,2,9,13)))        <- "Tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(35)))              <- "Npy+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(8,16)))            <- "Cell membrane projections"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(19)))              <- "Avp+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(25)))              <- "Oxt+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(0,12)))            <- "Astrocytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(18, 28)))          <- "Pars.Tuberalis"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(6)))               <- "Ependymocytes"

MEPV_cells$Level1 <- Idents(object = MEPV_cells)

plot1 = DimPlot(MEPV_cells, reduction = "humap") & theme(legend.text = element_text(size = 6)) & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = T)

################################################LEVEL2 annotation####################################
Idents(object = MEPV_cells) <- "seurat_clusters"

Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(4)))               <- "VLMC.1"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(21)))              <- "VLMC.2"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(29)))              <- "Dural fibroblasts"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(11)))              <- "Plvap Endothelial"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(1)))               <- "Endothelial cells"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(14, 24)))          <- "Pericytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(23)))              <- "VSMC"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(22, 31)))          <- "Immune cells"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(10, 30)))          <- "Microglia" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(17,32)))           <- "CAMs" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(15)))              <- "Progenitors"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(20)))              <- "Differentiating"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(33,26,3,5,34)))     <- "Mature"
#Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(27)))              <- "Lhb.Npy.Rax+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(7)))               <- "DMH tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(2)))               <- "VMH/dmARH tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(9)))               <- "vmARH tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(13)))              <- "ME tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(35)))              <- "Npy+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(8,16)))            <- "Cell membrane projections"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(19)))              <- "Avp+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(25)))              <- "Oxt+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(0,12)))            <- "Astrocytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(18)))              <- "Pars.Tuberalis"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(28)))              <- "Tight junction epithelial"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(6)))               <- "Ependymocytes"

MEPV_cells$Level2 <- Idents(object = MEPV_cells)
plot1 = DimPlot(MEPV_cells, reduction = "humap") & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = F) 

MEPV_cells
# An object of class Seurat 
# 49838 features across 47194 samples within 2 assays 
# Active assay: SCT (18914 features, 0 variable features)
# 1 other assay present: RNA
# 4 dimensional reductions calculated: pca, umap, harmony, humap


#No of cells in each cluster each condition
Idents (object = MEPV_cells) <-"Level2"
cell_num_v1 = table(Idents(MEPV_cells), MEPV_cells$Sample.name)
write.csv(cell_num_v1, "MEPV/output/Cell_number_v1.csv")

##Remove Lhb.Npy.Rax+ cells

Idents(object = MEPV_cells) <- "Level2"
MEPV_cells <-subset(MEPV_cells, idents=c("Lhb.Npy.Rax+"), invert = T)
MEPV_cells
# An object of class Seurat 
# 49838 features across 46922 samples within 2 assays 
# Active assay: SCT (18914 features, 0 variable features)
# 1 other assay present: RNA
# 4 dimensional reductions calculated: pca, umap, harmony, humap

#No of cells in each cluster each condition
Idents (object = MEPV_cells) <-"Level2"
cell_num_v2 = table(Idents(MEPV_cells), MEPV_cells$Sample.name)
write.csv(cell_num_v2, "MEPV/output/Cell_number_v2.csv")

##Rearrange Idents order 

MEPV_cells$Level2 <- factor(MEPV_cells$Level2,levels=c("Astrocytes", "Ependymocytes", "DMH tanycytes", "VMH/dmARH tanycytes", "vmARH tanycytes", 
                                                       "ME tanycytes", "Npy+", "Cell membrane projections", "Avp+", "Oxt+", 
                                                       "Microglia", "CAMs", "Immune cells", "Endothelial cells","Plvap Endothelial","Pericytes","VSMC",
                                                       "VLMC.1", "VLMC.2", "Dural fibroblasts", "Pars.Tuberalis", "Tight junction epithelial",
                                                       "Progenitors", "Differentiating", "Mature"))

###Assign colors to clusters####

col.pal <- list()
col.pal$celltype <- c("Astrocytes"="#A6CEE3", "DMH tanycytes"="#1F78B4", "VMH/dmARH tanycytes"="#B2DF8A", 
                      "vmARH tanycytes"="#33A02C", "ME tanycytes"="red",
                      "Npy+"="deeppink","Oxt+"="#FF7F00","Cell membrane projections"= "#CAB2D6", 
                      "Avp+"="#6A3D9A", "Immune cells"="#D95F02", "VLMC.1"="#abc4ff", 
                      "VLMC.2"="#B15928","Dural fibroblasts"="#CBD52E","Microglia"="#7570B3", 
                      "CAMs"="#E7298A", "Pars.Tuberalis"="#66A61E", "Tight junction epithelial"="#E6AB02",
                      "Plvap Endothelial"="#A6761D","Endothelial cells"="lightsalmon", "Pericytes"="yellow2","VSMC"="tomato1",
                      "Mature"="#0466c8", "Progenitors"="#38b000", "Differentiating"="ivory4", "Ependymocytes"="darkgoldenrod1")
col.pal <- list()
col.pal$celltype <- c("Astrocytes"="#829399", "DMH tanycytes"="#1F78B4", "VMH/dmARH tanycytes"="#B2DF8A", 
                      "vmARH tanycytes"="#33A02C", "ME tanycytes"="#FB9A99",
                      "Npy+"="#E31A1C","Oxt+"="#CAB2D6","Cell membrane projections"= "#6A3D9A", 
                      "Avp+"="#CBD52E", "Immune cells"="#B15928", "VLMC.1"="#1B9E77", 
                      "VLMC.2"="#D95F02","Dural fibroblasts"="#CBD52E","Microglia"="#7570B3", 
                      "CAMs"="#E7298A", "Pars.Tuberalis"="#66A61E", "Tight junction epithelial"="#E6AB02",
                      "Plvap Endothelial"="#FF7F00","Endothelial cells"="#abc4ff", "Pericytes"="yellow2","VSMC"="tomato1",
                      "Mature"="#0466c8", "Progenitors"="#38b000", "Differentiating"="ivory4", "Ependymocytes"="darkgoldenrod1")



pdf("MEPV/output/Fig1_Dimplot_harmony_labelled.pdf", width = 12.5, height = 8)
DimPlot(MEPV_cells, group.by = "seurat_clusters", cols = col.pal$celltype, reduction = "humap" , label.box = T) & NoLegend() & NoAxes()
dev.off()

pdf("MEPV/output/Fig1_Clustering_unlabelled.pdf", width = 6.1, height = 5)
DimPlot(MEPV_cells, group.by = "Level2", cols = (col.pal$celltype), reduction = "humap", pt.size = 0.7) & NoAxes()
dev.off()

pdf("MEPV/output/Fig1_Dimplot_harmony_clusters.pdf", width = 12.5, height = 8)
DimPlot(MEPV_cells, group.by = "seurat_clusters", reduction = "humap" , label = T) & NoAxes()
dev.off()

saveRDS(MEPV_cells, "MEPV/MEPV_cells_v3.rds")

##############Marker Plot################
Idents(MEPV_cells) <- "Level2"
Markers_ano <-FindAllMarkers(MEPV_cells, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.6)
Markers_ano2 <-FindAllMarkers(MEPV_cells, min.pct=0.25, only.pos = TRUE, logfc.threshold = 2)


write.csv(Markers_ano, "MEPV/output/Markers_ano2.csv")


Markers_ano2 %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
top5 <- top5[!duplicated(top5$gene),]

pdf("MEPV/output/Figure 1b_Dotplot_markers annotated_flip2.pdf", width = 9, height = 18)
DotPlot(object = MEPV_cells_1, 
        features = (top5$gene), 
        group.by = "Level2",
        assay = "SCT",
        scale = T,
        col.max = 2, 
        col.min = -2, cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression", 
                                 order = 1)) 
dev.off()

saveRDS(object = MEPV_cells, "MEPV/MEPV_cells_v2.rds")


#####Distribution of cell types across conditions#######

#reorder sample conditions
MEPV_cells$Sample.name <- factor(MEPV_cells$Sample.name,levels=c("Fem.Chow.Diest", "Fem.Chow.Proest", "Fem.Chow.Est", "Fem.HFDR", "Fem.HFDS", "Fem.HFDS+D", "Male.Chow", "Male.HFDR", "Male.HFDS"))

Idents(object = MEPV_cells) <- "Level2"
cells_by_cluster = as.matrix(table(Idents(MEPV_cells), MEPV_cells$Sample.name))

# 1. convert the data as a table
cells_by_cluster <- as.table(as.matrix(cells_by_cluster))
#Compute Chi-square residuals
chisq <- chisq.test(cells_by_cluster)
chisq

#Pearson's Chi-squared test
# data:  cells_by_cluster
# X-squared = 9007, df = 200, p-value < 2.2e-16

# Observed counts
chisq$observed

# Expected counts
round(chisq$expected,0)

round(chisq$residuals, 0)

library(corrplot)
pdf(file = "MEPV/output/Figure1_corrplot_cells_cond.pdf",   # The directory you want to save the file in
    width = 4.6, # The width of the plot in inches
    height = 6.2) # The height of the plot in inches
corrplot(chisq$residuals, is.cor = FALSE, tl.col = 'black', cl.cex = 0.8, tl.cex = 0.8,cl.ratio = 0.5,col = rev(brewer.pal(n=8, name="RdYlBu")))
dev.off()

#For a given cell, the size of the circle is proportional to the amount of the cell contribution.

# Contibution in percentage (%)
contrib <- 100*chisq$residuals^2/chisq$statistic
write.csv(contrib, "MEPV/output/contrib_cell_chi_sq.csv")
round(contrib, 3)

# Visualize the contribution
library(RColorBrewer)
corrplot(contrib, is.cor = FALSE, tl.col = 'black', cl.cex = 0.8, tl.cex = 0.8,cl.ratio = 0.5, col = (brewer.pal(n=8, name="YlOrRd")))

FeaturePlot(MEPV_cells, c("Tnfrsf11a"), reduction = "humap", order = TRUE, min.cutoff = "q50")
FeaturePlot(MEPV_cells, c("Tmem119"), reduction = "humap", order = TRUE, min.cutoff = "q50")
FeaturePlot(MEPV_cells, c("Tnfsf11"), reduction = "humap", order = TRUE)
FeaturePlot(MEPV_cells, c("Tppp3"), reduction = "humap", order = TRUE)

FeaturePlot(MEPV_cells, c("Olig1"), reduction = "humap", order = TRUE)

FeaturePlot(ventricuar_cells, c("IRF7"), reduction = "integrated", order = TRUE)

colnames(x = ventricuar_cells[["scvi"]]@cell.embeddings) <- paste0("scvi_", 1:2)

DimPlot(ventricuar_cells, group.by = "C4")

dittoBarPlot(MEPV_cells, "Level2", group.by = "Sample.name")

FeaturePlot(MEPV_cells, c("Irf7"), reduction = "humap", order = TRUE ) & scale_colour_viridis(option = "C", direction = -1, na.value = "grey50")

library(dittoSeq)
pdf("MEPV/output/Fig1_DittoMEPV_labelled.pdf", width = 12.5, height = 8)
dittoDimPlot(MEPV_cells, var = "Level2",reduction.use = 'humap', opacity = 0.5)
dev.off()


#proportion of OSNs by PCW
levels(Tanycytes_sub) <- c("Fem.Chow.Proest", "Fem.Chow.Est", "Fem.Chow.Diest", "Fem.HFDR", "Fem.HFDS", "Fem.HFDS+D", "Male.Chow", "Male.HFDR", "Male.HFDS")


Idents(Tanycytes_sub) <-  "Tany.cell.state"
Tany8 <-subset(Tanycytes_sub, idents=c("Tany.8"))


saveRDS(object = MEPV_cells, "D:/PUBLICATIONS AND REVIEWS/Manuscript_single_cell_2024/Manu_final_submit/Raw data/Script/MEPV_cells.rds")
