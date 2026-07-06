############################################################
## Final Human tanycyte annotation using:
## 1. Human and mouse marker overlap
## 2. Jaccard similarity
## 3. scANVI transferred labels
## 4. Annotation table and visualizations
############################################################

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(proxy)
  library(scales)
  })

############################################################
## 0. Paths and settings
############################################################

base_dir <- "D:/PROJECTS/snRNseq_integration/human_ME_full"
out_dir  <- file.path(base_dir, "Part3_cross_species", "after_scvi")
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

human_rds <- file.path(out_dir, "Tanycytes_human_with_scanvi.rds")
mouse_rds <- file.path(base_dir, "Part1_Mouse_tany", "Tanycytes_sub.rds")

scanvi_label_column <- "scanvi_final_label"
human_top_n_markers <- 250
mouse_top_n_markers <- 200



############################################################
## 1. Load objects
############################################################

human <- readRDS(human_rds)
mouse <- readRDS(mouse_rds)

DefaultAssay(human) <- "SCT"
DefaultAssay(mouse) <- "SCT"

#stopifnot(scanvi_label_column %in% colnames(human@meta.data))

############################################################
## 2. Helper functions
############################################################

get_top_markers <- function(obj, species, top_n, ident_col = NULL) {
  if (!is.null(ident_col)) {
    Idents(obj) <- ident_col
  }
  
  FindAllMarkers(
    obj,
    min.pct = 0.25,
    only.pos = TRUE,
    logfc.threshold = 0.1
  ) %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
    mutate(
      gene = toupper(gene),
      species = species,
      cluster_name = paste0(species, "_", cluster)
    ) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup()
}

make_marker_sets <- function(marker_df) {
  marker_df %>%
    group_by(cluster_name) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
}

make_overlap_table <- function(human_marker_sets, mouse_marker_sets) {
  overlap_table <- expand.grid(
    human_cluster = human_marker_sets$cluster_name,
    mouse_cluster = mouse_marker_sets$cluster_name,
    stringsAsFactors = FALSE
  )
  
  overlap_table$overlap_n <- NA_integer_
  overlap_table$jaccard <- NA_real_
  overlap_table$shared_gene_string <- NA_character_
  
  for (i in seq_len(nrow(overlap_table))) {
    h <- overlap_table$human_cluster[i]
    m <- overlap_table$mouse_cluster[i]
    
    h_genes <- human_marker_sets$genes[human_marker_sets$cluster_name == h][[1]]
    m_genes <- mouse_marker_sets$genes[mouse_marker_sets$cluster_name == m][[1]]
    shared <- intersect(h_genes, m_genes)
    
    overlap_table$overlap_n[i] <- length(shared)
    overlap_table$jaccard[i] <- length(shared) / length(union(h_genes, m_genes))
    overlap_table$shared_gene_string[i] <- paste(shared, collapse = ", ")
  }
  
  overlap_table %>%
    arrange(human_cluster, desc(overlap_n), desc(jaccard))
}

make_jaccard_matrix <- function(overlap_table) {
  overlap_table %>%
    select(human_cluster, mouse_cluster, jaccard) %>%
    pivot_wider(names_from = mouse_cluster, values_from = jaccard) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("human_cluster") %>%
    as.matrix()
}

make_overlap_matrix <- function(overlap_table) {
  overlap_table %>%
    select(human_cluster, mouse_cluster, overlap_n) %>%
    pivot_wider(names_from = mouse_cluster, values_from = overlap_n) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("human_cluster") %>%
    as.matrix()
}
############################################################
## 3. Marker detection and marker sets
############################################################
###Analysis with resolution 1.2
# Idents(human) <- "SCT_snn_res.1.2"
# human$seurat_clusters <- human$SCT_snn_res.1.2
Idents(human) <- "SCT_snn_res.1.2"
human$seurat_clusters <- as.character(human$SCT_snn_res.1.2)
Idents(human) <- "seurat_clusters"

human_markers <- get_top_markers(
  obj = human,
  species = "Human",
  top_n = human_top_n_markers,
  ident_col = "seurat_clusters"
)

mouse_markers <- get_top_markers(
  obj = mouse,
  species = "Mouse",
  top_n = mouse_top_n_markers
)

human_marker_sets <- make_marker_sets(human_markers)
mouse_marker_sets <- make_marker_sets(mouse_markers)

############################################################
## 4. Human-mouse marker overlap and Jaccard similarity
############################################################

all_overlap_table <- make_overlap_table(human_marker_sets, mouse_marker_sets)
jaccard_mat <- make_jaccard_matrix(all_overlap_table)
overlap_mat <- make_overlap_matrix(all_overlap_table)

best_marker_match <- all_overlap_table %>%
  group_by(human_cluster) %>%
  arrange(desc(jaccard), desc(overlap_n), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    human_cluster,
    best_mouse_marker_match = mouse_cluster,
    jaccard_score = jaccard,
    overlap_genes = overlap_n,
    shared_gene_string
  )

############################################################
## 5. scANVI transferred-label summary
############################################################

scanvi_summary <- human@meta.data %>%
  mutate(
    human_cluster = paste0("Human_", seurat_clusters),
    scanvi_label = .data[[scanvi_label_column]]
  ) %>%
  count(human_cluster, scanvi_label, name = "scanvi_n") %>%
  group_by(human_cluster) %>%
  mutate(scanvi_fraction = scanvi_n / sum(scanvi_n)) %>%
  arrange(human_cluster, desc(scanvi_fraction)) %>%
  ungroup()

best_scanvi_per_human <- scanvi_summary %>%
  group_by(human_cluster) %>%
  slice_max(scanvi_fraction, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    dominant_scanvi_label = scanvi_label,
    dominant_scanvi_n = scanvi_n,
    dominant_scanvi_fraction = scanvi_fraction
  )

scanvi_all_labels <- scanvi_summary %>%
  mutate(
    scanvi_label_summary = paste0(
      scanvi_label,
      " (", scanvi_n, ", ", round(scanvi_fraction * 100, 1), "%)"
    )
  ) %>%
  group_by(human_cluster) %>%
  summarise(
    all_scanvi_labels = paste(scanvi_label_summary, collapse = "; "),
    .groups = "drop"
  )
############################################################
## 6. Combined annotation/support tables
############################################################

all_overlap_table_scanvi <- all_overlap_table %>%
  left_join(best_scanvi_per_human, by = "human_cluster") %>%
  left_join(scanvi_all_labels, by = "human_cluster") %>%
  mutate(
    mouse_cluster_clean = gsub("^Mouse_", "", mouse_cluster),
    scanvi_marker_agreement = mouse_cluster_clean == dominant_scanvi_label
  ) %>%
  arrange(human_cluster, desc(overlap_n), desc(jaccard))

annotation_table <- best_marker_match %>%
  left_join(best_scanvi_per_human, by = "human_cluster") %>%
  mutate(
    marker_mouse_label = gsub("^Mouse_", "", best_mouse_marker_match),
    marker_scanvi_agreement = marker_mouse_label == dominant_scanvi_label,
    suggested_annotation = ifelse(
      marker_scanvi_agreement,
      paste0(dominant_scanvi_label, "-like"),
      paste0(marker_mouse_label, "-like / check manually")
    )
  )

write.csv(
  annotation_table,
  file.path(out_dir, "human_mouse_tanycyte_annotation_table.csv"),
  row.names = FALSE
)

write.csv(
  all_overlap_table_scanvi,
  file.path(out_dir, "human_mouse_tanycyte_all_overlap_with_scanvi.csv"),
  row.names = FALSE
)

write.csv(
  scanvi_summary,
  file.path(out_dir, "human_scanvi_label_summary_by_cluster.csv"),
  row.names = FALSE
)


## Use all pairwise human-mouse comparisons
bubble_df <- all_overlap_table_scanvi %>%
  mutate(
    mouse_cluster = gsub("^Mouse_", "", mouse_cluster),
    scanvi_agreement = ifelse(
      scanvi_marker_agreement,
      "Agreement",
      "Discordant"
    )
  )

## Add final/best marker match flag
best_pairs <- bubble_df %>%
  group_by(human_cluster) %>%
  arrange(desc(jaccard), desc(overlap_n), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(human_cluster, mouse_cluster) %>%
  mutate(final_best_marker_match = TRUE)

bubble_df <- bubble_df %>%
  left_join(best_pairs, by = c("human_cluster", "mouse_cluster")) %>%
  mutate(
    final_best_marker_match = ifelse(
      is.na(final_best_marker_match),
      FALSE,
      final_best_marker_match
    )
  )

## order mouse clusters 
bubble_df <- bubble_df %>%
  mutate(
    mouse_cluster_clean = factor(
      mouse_cluster_clean,
      levels = paste0("mT.", 1:8)
    )
  )

## order human clusters by their best mouse marker match
human_order <- bubble_df %>%
  filter(final_best_marker_match) %>%
  arrange(mouse_cluster, human_cluster) %>%
  pull(human_cluster)

bubble_df <- bubble_df %>%
  mutate(
    human_cluster = factor(
      human_cluster,
      levels = rev(unique(human_order))
    )
  )

bubble_df <- bubble_df %>%
  mutate(
    mouse_label_clean = gsub("^Mouse_", "", mouse_cluster),
    scanvi_predicted_match = mouse_label_clean == dominant_scanvi_label
  )

p_all_marker_scanvi <- ggplot(
  bubble_df,
  aes(x = mouse_cluster, y = human_cluster)
) +
  ## Base bubbles: all marker comparisons
  geom_point(
    aes(
      size = overlap_n,
      fill = jaccard
    ),
    shape = 21,
    color = "grey75",
    stroke = 0.7,
    alpha = 0.95
  ) +
  
  ## Outer blue ring: dominant scANVI-predicted mouse label
  geom_point(
    data = bubble_df %>% filter(scanvi_predicted_match),
    aes(
      size = overlap_n,
      fill = jaccard
    ),
    shape = 21,
    color = "#0072B2",
    stroke = 2.0,
    alpha = 0.95
  ) +
  
  ## Inner black ring: best marker match
  geom_point(
    data = bubble_df %>% filter(final_best_marker_match),
    aes(
      size = overlap_n,
      fill = jaccard
    ),
    shape = 21,
    color = "black",
    stroke = 0.8,
    alpha = 0.95
  ) +
  
  scale_size_area(
    max_size = 10,
    name = "Shared marker genes"
  ) +
  scale_fill_gradientn(
    colours = c("white", "gold", "darkred"),
    name = "Jaccard similarity"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.box = "vertical"
  ) +
  labs(
    title = "Human-mouse tanycyte marker similarity with scANVI support",
    subtitle = "Size = shared markers; fill = Jaccard; blue ring = dominant scANVI label",
    x = "Mouse tanycyte subtype",
    y = "Human cluster"
  )

p_all_marker_scanvi

ggsave(
  filename = file.path(out_dir, "human_mouse_marker_scanvi_bubble_final.pdf"),
  plot = p_all_marker_scanvi,
  width = 5.7,
  height = 6.0,
  units = "in"
)

############################################################
## 9. Final human subtype annotations
############################################################

##First Dimplot with seurat clusters

human_cols_seurat <- setNames(
  c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "#A6CEE3", "#FB9A99", "#CAB2D6", "#B2DF8A",
    "#FDBF6F", "#B15928"
  ),
  paste0("Human_", 0:13)
)

human$human_cluster <- paste0("Human_", human$seurat_clusters)

human$human_cluster <- factor(
  human$human_cluster,
  levels = paste0("Human_", 0:13)
)

Idents(human) <- "human_cluster"
p_human_clusters_umap <- DimPlot(
  human,
  reduction = "umap",
  group.by = "human_cluster",
  cols = human_cols_seurat,
  label = FALSE,
  order = TRUE,
  pt.size = 0.1
) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("Human tanycyte Seurat clusters")

p_human_clusters_umap

ggsave(
  filename = file.path(out_dir, "human_seurat_clusters_umap.pdf"),
  plot = p_human_clusters_umap,
  width = 6.3,
  height = 5.1,
  units = "in"
)


###Final annotations
Idents(human) <- "seurat_clusters"
Idents(human, cells = WhichCells(human, idents = c(3,6,5,9))) <- "hT.1"
Idents(human, cells = WhichCells(human, idents = c(8))) <- "hT.3"
Idents(human, cells = WhichCells(human, idents = c(12))) <- "hT.4"
Idents(human, cells = WhichCells(human, idents = c(2,13))) <- "hT.5-6"
Idents(human, cells = WhichCells(human, idents = c(0,1,10))) <- "hT.7"
# Idents(human, cells = WhichCells(human, idents = c(11))) <- "hT.8_1"
Idents(human, cells = WhichCells(human, idents = c(4,7,11))) <- "hT.8"

human$Tany_grouped <- factor(
  Idents(human),
  levels = c("hT.1","hT.3", "hT.4","hT.5-6", "hT.7", "hT.8")
)

Idents(human) <- "Tany_grouped"

human_cols_ano <- c(
  "hT.1"   = "#E69F00",
  #"hT.2"   = "#56B4E9",   # blue
  "hT.3"   = "#009E73",
  "hT.4"   = "#F0E442",
  "hT.5-6"   = "#D55E00",
  "hT.7"   = "#F781BF",
  #"hT.8_1" = "#984EA3",
  "hT.8" = "#9E9E9E"
)

p_human <- DimPlot(
  human,
  group.by = "Tany_grouped",
  cols = human_cols_ano,
  label = FALSE,
  order = TRUE,
  pt.size = 0.1
) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank()
  ) 

print(p_human)
ggsave(
  filename = file.path(out_dir, "Dim_Human_tanycyte_subtypes.pdf"),
  plot = p_human,
  width =5.6,
  height = 4.5,
  units = "in"
)

saveRDS(
  human,
  file.path(out_dir, "Tanycytes_human_with_scanvi_grouped_ano.rds")
)

########Final dendrogram based on new annotations########

############################################################
## 1. Set final human annotations
############################################################

Idents(human) <- "Tany_grouped"

DefaultAssay(human) <- "SCT"
DefaultAssay(mouse) <- "SCT"

############################################################
## 2. Set mouse annotations
############################################################
## Replace this if your mouse subtype column has a different name

Idents(mouse) <- "Tany.rename"

############################################################
## 3. Find markers for final human tanycyte groups
############################################################

human_group_markers <- FindAllMarkers(
  human,
  min.pct = 0.25,
  only.pos = TRUE,
  logfc.threshold = 0.1
) %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  mutate(
    gene = toupper(gene),
    cluster_name = as.character(cluster),
    dendro_name = paste0("Human_", cluster_name)
  )

top_human_markers <- human_group_markers %>%
  group_by(dendro_name) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()

############################################################
## 4. Find markers for mouse tanycyte groups
############################################################

mouse_group_markers <- FindAllMarkers(
  mouse,
  min.pct = 0.25,
  only.pos = TRUE,
  logfc.threshold = 0.1
) %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  mutate(
    gene = toupper(gene),
    cluster_name = as.character(cluster),
    dendro_name = paste0("Mouse_", cluster_name)
  )

top_mouse_markers <- mouse_group_markers %>%
  group_by(dendro_name) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()

############################################################
## 5. Build marker gene sets
############################################################

human_marker_sets_final <- top_human_markers %>%
  group_by(dendro_name) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  )

mouse_marker_sets_final <- top_mouse_markers %>%
  group_by(dendro_name) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  )

all_marker_sets_final <- bind_rows(
  human_marker_sets_final,
  mouse_marker_sets_final
)

all_genes_final <- sort(unique(unlist(all_marker_sets_final$genes)))

############################################################
## 6. Binary cluster × gene matrix
############################################################

marker_binary_final <- matrix(
  0,
  nrow = nrow(all_marker_sets_final),
  ncol = length(all_genes_final),
  dimnames = list(all_marker_sets_final$dendro_name, all_genes_final)
)

for (i in seq_len(nrow(all_marker_sets_final))) {
  marker_binary_final[
    all_marker_sets_final$dendro_name[i],
    all_marker_sets_final$genes[[i]]
  ] <- 1
}

############################################################
## 7. Jaccard dendrogram
############################################################

cluster_dist_final <- proxy::dist(
  marker_binary_final,
  method = "Jaccard"
)

cluster_tree_final <- hclust(
  cluster_dist_final,
  method = "average"
)

plot(
  cluster_tree_final,
  main = "Human-mouse tanycyte marker dendrogram",
  xlab = "",
  sub = "",
  cex = 0.9
)

############################################################
## 8. Save dendrogram as PDF
############################################################

pdf(
  file = file.path(out_dir, "human_mouse_final_tanycyte_marker_dendrogram.pdf"),
  width = 6.8,
  height = 4.6
)

plot(
  cluster_tree_final,
  main = "Human-mouse tanycyte marker dendrogram",
  xlab = "",
  sub = "",
  cex = 0.9
)

dev.off()

######Human tanycyte marker genes plot#########
human_genes_to_plot2 = c("GPR37L1", "NTSR2", "RFX4", "SEMA6D","SLC1A2", "SHISA9", "RGS20", #hT.1 and mT.1
                        "RSPO3","TGFB2","FLT1","NR2F2", #mT2 unique; hT2 not present in humans
                        "ADGB", "VWA5B1", "DNAH11", #present in mT3 and mT4  
                        "TPPP3", "CTXN1", "HIPK1", #present in both mT3 and mT4 but enriched in hT4 
                        "VCAN","CRYM","EPHB1","RORB",
                        "FGF10", "TOX2", "MOB3B","FRZB","ADAMTSL1","GRIA3","DIO2","GPC3","DEPTOR", 
                        "COL25A1", "MEST", "A2M","B2M", #ME T7
                        "CD44","IFI44", "IFITM2", "IFITM3","IFI27","GPX1",  #immune - T8
                        "CITED1", "HLA-C", "H3-3A","CD59","CD74") # "CXCL14", "IFI27","HLA-DRB1") 


human$Tany_grouped <- factor(
  human$Tany_grouped,
  levels = rev(c(
    "hT.1",
    "hT.3",
    "hT.4",
    "hT.5-6",
    "hT.7",
    "hT.8"
  ))
)

Idents(human) <- "Tany_grouped"
dot_human = DotPlot(object = human, ,
        features = human_genes_to_plot2,
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
 
dot_human
ggsave(
  filename = file.path(plot_dir, "Dotplot_human_markers.pdf"),
  plot = dot_human,
  width = 13.7,
  height = 3.2,
  units = "in"
)
