####################################
# GO Biological Process enrichment
####################################

library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

set.seed(1234)

out_dir <- "human_ME_full/Part3_cross_species/after_scvi"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------
# 1. GO enrichment function
#-----------------------------

run_go_bp <- function(marker_df, orgdb, top_n = 500) {
  
  gene_lists <- marker_df %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0, !is.na(gene)) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = top_n, with_ties = FALSE) %>%
    summarise(
      genes = list(unique(gene)),
      .groups = "drop"
    )
  
  gene_lists <- setNames(gene_lists$genes, gene_lists$cluster)
  
  go_results <- lapply(
    gene_lists,
    function(genes) {
      enrichGO(
        gene = genes,
        OrgDb = orgdb,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH"
      )
    }
  )
  
  bind_rows(
    lapply(names(go_results), function(cl) {
      res <- as.data.frame(go_results[[cl]])
      
      if (nrow(res) == 0) {
        return(NULL)
      }
      
      res %>% mutate(cluster = cl)
    })
  )
}

#-----------------------------
# 2. Run enrichment
#-----------------------------

Idents(human) <- "Tany_grouped"
Idents(mouse) <- "Tany.cell.state"

human_go_table <- run_go_bp(
  marker_df = human_group_markers,
  orgdb = org.Hs.eg.db,
  top_n = 500
)

write.csv(
  human_go_table,
  file.path(out_dir, "human_GO_BP_all_results_500.csv"),
  row.names = FALSE
)

mouse_group_markers <- FindAllMarkers(
  mouse,
  min.pct = 0.25,
  only.pos = TRUE,
  logfc.threshold = 0.1
)

mouse_go_table <- run_go_bp(
  marker_df = mouse_group_markers,
  orgdb = org.Mm.eg.db,
  top_n = 500
)

write.csv(
  mouse_go_table,
  file.path(out_dir, "mouse_GO_BP_all_results_500.csv"),
  row.names = FALSE
)

#-----------------------------
# 3. Read manually selected GO terms
#-----------------------------

go_df <- read_csv(
  file.path(out_dir, "Human_mouse_GO_BP_500_selected.csv"),
  show_col_types = FALSE
)

cluster_order <- c(
  paste0("mT.", 1:8),
  "hT.1", "hT.3", "hT.4", "hT.5-6", "hT.7", "hT.8"
)

cluster_colors <- c(
  "mT.1" = "#E69F00",
  "mT.2" = "#56B4E9",
  "mT.3" = "#009E73",
  "mT.4" = "#F0E442",
  "mT.5" = "#0072B2",
  "mT.6" = "#D55E00",
  "mT.7" = "#F781BF",
  "mT.8" = "#9E9E9E",
  "hT.1" = "#E69F00",
  "hT.3" = "#009E73",
  "hT.4" = "#F0E442",
  "hT.5-6" = "#D55E00",
  "hT.7" = "#F781BF",
  "hT.8" = "#9E9E9E"
)

#-----------------------------
# 4. Build GO term x cluster matrix
#-----------------------------

heat_df <- go_df %>%
  group_by(Description, cluster) %>%
  summarise(
    zScore = max(zScore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = zScore,
    values_fill = 0
  )

heat_mat <- as.matrix(heat_df[, -1])
rownames(heat_mat) <- heat_df$Description

heat_mat <- heat_mat[
  ,
  intersect(cluster_order, colnames(heat_mat)),
  drop = FALSE
]

#-----------------------------
# 5. Scale GO enrichment within species
#-----------------------------

scale_0_1 <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(rep(0, length(x)))
  }
  
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

mouse_cols <- grep("^mT\\.", colnames(heat_mat), value = TRUE)
human_cols <- grep("^hT\\.", colnames(heat_mat), value = TRUE)

mouse_scaled <- t(apply(heat_mat[, mouse_cols, drop = FALSE], 1, scale_0_1))
human_scaled <- t(apply(heat_mat[, human_cols, drop = FALSE], 1, scale_0_1))

scaled_mat <- cbind(mouse_scaled, human_scaled)

#-----------------------------
# 6. Prepare long-format plot table
#-----------------------------

plot_df <- as.data.frame(scaled_mat) %>%
  rownames_to_column("Description") %>%
  pivot_longer(
    cols = -Description,
    names_to = "cluster",
    values_to = "scaled_score"
  ) %>%
  mutate(
    Species = ifelse(grepl("^mT\\.", cluster), "Mouse", "Human"),
    Species = factor(Species, levels = c("Mouse", "Human")),
    cluster = factor(cluster, levels = cluster_order)
  )

term_order <- plot_df %>%
  group_by(Description) %>%
  slice_max(scaled_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(scaled_score)) %>%
  pull(Description)

plot_df <- plot_df %>%
  mutate(
    Description_short = str_wrap(Description, width = 45),
    Description_short = factor(
      Description_short,
      levels = rev(str_wrap(term_order, width = 45))
    )
  )

range_df <- plot_df %>%
  group_by(Description_short, Species) %>%
  summarise(
    xmin = min(scaled_score, na.rm = TRUE),
    xmax = max(scaled_score, na.rm = TRUE),
    .groups = "drop"
  )

label_df <- plot_df %>%
  group_by(Description_short, Species) %>%
  slice_max(scaled_score, n = 1, with_ties = FALSE) %>%
  ungroup()

#-----------------------------
# 7. Dumbbell-style GO profile plot
#-----------------------------

p_go <- ggplot() +
  geom_segment(
    data = range_df,
    aes(
      x = xmin,
      xend = xmax,
      y = Description_short,
      yend = Description_short
    ),
    color = "grey90",
    linewidth = 0.45
  ) +
  geom_point(
    data = plot_df,
    aes(
      x = scaled_score,
      y = Description_short,
      color = cluster
    ),
    size = 2.3,
    alpha = 0.9
  ) +
  geom_text(
    data = label_df,
    aes(
      x = scaled_score,
      y = Description_short,
      label = cluster
    ),
    hjust = -0.25,
    size = 2.3
  ) +
  facet_grid(
    ~ Species,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_manual(
    values = cluster_colors,
    breaks = cluster_order,
    name = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Relative enrichment within species",
    y = "",
    title = "Cluster-level GO Biological Process profiles"
  )

ggsave(
  filename = file.path(plot_dir, "GO_BP_species_minmax_dumbbell_profile.pdf"),
  plot = p_go,
  width = 6.5,
  height = 5.3
)

p_go







































####################################
######Pathway analysis##############
####################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

Idents(human) <- "Tany_grouped"

human_top_enrich_genes <- human_group_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 500)

human_enrich_gene_lists <- human_top_enrich_genes %>%
  filter(!is.na(gene)) %>%
  mutate(gene = toupper(gene)) %>%
  group_by(cluster) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  )

human_enrich_gene_lists <- setNames(
  human_enrich_gene_lists$genes,
  human_enrich_gene_lists$cluster
)

human_go <- lapply(
  human_enrich_gene_lists,
  function(x)
    enrichGO(
      gene = x,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH"
    )
)

##save
human_go_table <- bind_rows(
  lapply(names(human_go), function(cl){
    
    res <- human_go[[cl]]
    
    if(is.null(res) || nrow(as.data.frame(res)) == 0)
      return(NULL)
    
    as.data.frame(res) %>%
      mutate(cluster = cl)
    
  })
)

write.csv(
  human_go_table,
  file.path(out_dir, "human_GO_BP_all_results_500.csv"),
  row.names = FALSE
)


###Mouse GOBPs######
mouse_group_markers_scase <- FindAllMarkers(
  mouse,
  min.pct = 0.25,
  only.pos = TRUE,
  logfc.threshold = 0.1
)


mouse_top_enrich_genes <- mouse_group_markers_scase %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 500)

mouse_enrich_gene_lists <- mouse_top_enrich_genes %>%
  filter(!is.na(gene)) %>%
  group_by(cluster) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  )

mouse_enrich_gene_lists <- setNames(
  mouse_enrich_gene_lists$genes,
  mouse_enrich_gene_lists$cluster
)

mouse_go <- lapply(
  mouse_enrich_gene_lists,
  function(x)
    enrichGO(
      gene = x,
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH"
    )
)

mouse_go_table <- bind_rows(
  lapply(names(mouse_go), function(cl){
    
    res <- mouse_go[[cl]]
    
    if(is.null(res) || nrow(as.data.frame(res)) == 0)
      return(NULL)
    
    as.data.frame(res) %>%
      mutate(cluster = cl)
    
  })
)

write.csv(
  mouse_go_table,
  file.path(out_dir, "mouse_GO_BP_all_results_500.csv"),
  row.names = FALSE
)

library(dplyr)
library(tidyr)
library(pheatmap)

############################################################
## 1. Prepare human and mouse GO tables
############################################################

human_go_plot <- human_go_table %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    species_cluster = paste0("Human_", cluster),
    score = -log10(p.adjust)
  )

mouse_go_plot <- mouse_go_table %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    species_cluster = paste0("Mouse_", cluster),
    score = -log10(p.adjust)
  )

############################################################
## 2. Find common GO-BP terms
############################################################

common_terms <- intersect(
  unique(human_go_plot$Description),
  unique(mouse_go_plot$Description)
)

length(common_terms)

common_go_table <- bind_rows(
  human_go_plot,
  mouse_go_plot
) %>%
  filter(Description %in% common_terms)

top_common_terms <- common_go_table %>%
  group_by(Description) %>%
  summarise(
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_score)) %>%
  slice_head(n = 50) %>%
  pull(Description)

go_heat_df <- common_go_table %>%
  filter(Description %in% top_common_terms) %>%
  dplyr::select(Description, species_cluster, score) %>%
  group_by(Description, species_cluster) %>%
  summarise(
    score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = species_cluster,
    values_from = score,
    values_fill = 0
  )

go_heat_mat <- as.matrix(go_heat_df[, -1])
rownames(go_heat_mat) <- go_heat_df$Description

human_order <- paste0(
  "Human_",
  c("hT.1", "hT.3", "hT.4", "hT.5-6", "hT.7", "hT.8")
)

mouse_order <- paste0(
  "Mouse_",
  paste0("mT.", 1:8)
)

Human_mouse_GO_BP_500_selected <- read_csv("human_ME_full/Part3_cross_species/after_scvi/Human_mouse_GO_BP_500_selected.csv")
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)
library(ggrepel)

go_df <- Human_Mouse_GO_BP_common_500_selected

cluster_order <- c(
  paste0("mT.", 1:8),
  "hT.1", "hT.3", "hT.4", "hT.5-6", "hT.7", "hT.8"
)

#########################################################
## Build GO term x cluster matrix using original zScore
#########################################################

heat_df <- go_df %>%
  group_by(Description, cluster) %>%
  summarise(
    zScore = max(zScore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = zScore,
    values_fill = 0
  )

heat_mat <- as.matrix(heat_df[, -1])
rownames(heat_mat) <- heat_df$Description

heat_mat <- heat_mat[
  ,
  intersect(cluster_order, colnames(heat_mat)),
  drop = FALSE
]

#########################################################
## Split mouse and human
#########################################################

mouse_cols <- grep("^mT\\.", colnames(heat_mat), value = TRUE)
human_cols <- grep("^hT\\.", colnames(heat_mat), value = TRUE)

mouse_mat <- heat_mat[, mouse_cols, drop = FALSE]
human_mat <- heat_mat[, human_cols, drop = FALSE]

#########################################################
## Min-max scale each GO term separately within species
#########################################################

scale_0_1 <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    rep(0, length(x))
  } else {
    (x - min(x, na.rm = TRUE)) /
      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
}

mouse_scaled <- t(apply(mouse_mat, 1, scale_0_1))
human_scaled <- t(apply(human_mat, 1, scale_0_1))

colnames(mouse_scaled) <- colnames(mouse_mat)
colnames(human_scaled) <- colnames(human_mat)
rownames(mouse_scaled) <- rownames(mouse_mat)
rownames(human_scaled) <- rownames(human_mat)

scaled_mat <- cbind(mouse_scaled, human_scaled)

#########################################################
## Long format
#########################################################

plot_df <- as.data.frame(scaled_mat) %>%
  rownames_to_column("Description") %>%
  pivot_longer(
    cols = -Description,
    names_to = "cluster",
    values_to = "scaled_score"
  ) %>%
  mutate(
    Species = ifelse(grepl("^mT\\.", cluster), "Mouse", "Human"),
    Species = factor(Species, levels = c("Mouse", "Human")),
    cluster = factor(cluster, levels = cluster_order)
  )

#########################################################
## Row order by dominant cluster
#########################################################

term_order <- plot_df %>%
  group_by(Description) %>%
  slice_max(order_by = scaled_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(scaled_score)) %>%
  pull(Description)

plot_df <- plot_df %>%
  mutate(
    Description_short = str_wrap(Description, width = 45),
    Description_short = factor(
      Description_short,
      levels = rev(str_wrap(term_order, width = 45))
    )
  )

#########################################################
## Get min/max and dominant cluster per GO term/species
#########################################################

range_df <- plot_df %>%
  group_by(Description_short, Species) %>%
  summarise(
    xmin = min(scaled_score, na.rm = TRUE),
    xmax = max(scaled_score, na.rm = TRUE),
    .groups = "drop"
  )

label_df <- plot_df %>%
  group_by(Description_short, Species) %>%
  slice_max(order_by = scaled_score, n = 1, with_ties = FALSE) %>%
  ungroup()

cluster_colors <- c(
  "mT.1" = "#E69F00",
  "mT.2" = "#56B4E9",
  "mT.3" = "#009E73",
  "mT.4" = "#F0E442",
  "mT.5" = "#0072B2",   # darker blue / dittoSeq-like
  "mT.6" = "#D55E00",
  "mT.7" = "#F781BF",
  "mT.8" = "#9E9E9E",
  
  "hT.1" = "#E69F00",
  "hT.3" = "#009E73",
  "hT.4" = "#F0E442",
  "hT.5-6" = "#D55E00",
  "hT.7" = "#F781BF",
  "hT.8" = "#9E9E9E"
)

scale_color_manual(values = cluster_colors)

#########################################################
## Dumbbell-style profile plot
#########################################################

p <- ggplot() +
  geom_segment(
    data = range_df,
    aes(
      x = xmin,
      xend = xmax,
      y = Description_short,
      yend = Description_short
    ),
    color = "grey90",
    linewidth = 0.45
  ) +
  geom_point(
    data = plot_df,
    aes(
      x = scaled_score,
      y = Description_short,
      color = cluster
    ),
    size = 2.3,
    alpha = 0.9
  ) +
  geom_text(
    data = label_df,
    aes(
      x = scaled_score,
      y = Description_short,
      label = cluster
    ),
    hjust = -0.25,
    size = 2.3
  ) +
  facet_grid(
    ~ Species,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_manual(
    values = cluster_colors,
    breaks = cluster_order,
    name = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Relative enrichment within species",
    y = "",
    title = "Cluster-level GO Biological Process profiles"
  )

p

ggsave(
  file.path(plot_dir, "GO_BP_species_minmax_dumbbell_profile.pdf"),
  plot = p,
  width = 6.5,
  height = 5.3
)

