# UMAP Plots --------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Generates UMAP visualizations using lightweight data extracted
#              from the processed Seurat object. Includes UMAPs colored by
#              sample, sex, cluster, percent MT, and annotated cell type.
#              Cell type annotations are first-pass based on top marker genes
#              and known ileum biology — refine as needed.
#              Style: no titles, no captions, minimal arrow axes.
# Input:       data/processed/umap_embeddings.rds
#              data/processed/metadata_processed.rds
#              data/processed/elbow_df.rds
# Output:      plots/pca_elbow_plot.pdf
#              plots/umap_by_sample.pdf
#              plots/umap_by_sex.pdf
#              plots/umap_clusters.pdf
#              plots/umap_pct_mt.pdf
#              plots/umap_annotated.pdf
#              plots/umap_annotated_labeled.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Default clustering resolution — must match normalize_cluster.R
default_resolution <- 0.5

## Cluster column name in metadata
cluster_col <- paste0("SCT_snn_res.", default_resolution)

## First-pass cluster annotations
cluster_annotations <- c(
  "0"  = "Enterocyte",
  "1"  = "Enteric Neuron",
  "2"  = "Unknown",
  "3"  = "T cell",
  "4"  = "Enterocyte",
  "5"  = "Fibroblast",
  "6"  = "Enteric Neuron",
  "7"  = "Immune",
  "8"  = "Enterocyte",
  "9"  = "Goblet",
  "10" = "T cell",
  "11" = "Smooth Muscle",
  "12" = "Endothelial",
  "13" = "Lymphatic Endothelial",
  "14" = "Enteric Neuron",
  "15" = "Fibroblast",
  "16" = "Enteric Neuron",
  "17" = "Enteric Neuron",
  "18" = "Enteric Neuron",
  "19" = "Enteric Glia",
  "20" = "Enteric Neuron",
  "21" = "Enteric Neuron",
  "22" = "Enteric Neuron",
  "23" = "Enteric Glia",
  "24" = "Smooth Muscle",
  "25" = "Smooth Muscle",
  "26" = "Smooth Muscle",
  "27" = "Unknown",
  "28" = "Smooth Muscle"
)

## Cell type color palette
celltype_colors <- c(
  "Enterocyte"            = "#4E9AF1",
  "Goblet"                = "#6DC8A0",
  "T cell"                = "#06D6A0",
  "Immune"                = "#2EC4B6",
  "Enteric Neuron"        = "#C77DFF",
  "Enteric Glia"          = "#F15BB5",
  "Smooth Muscle"         = "#FF6B6B",
  "Endothelial"           = "#FFD166",
  "Lymphatic Endothelial" = "#F4A261",
  "Fibroblast"            = "#118AB2",
  "Unknown"               = "grey70"
)

## Arrow length in UMAP coordinate units
arrow_length <- 2


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(here)


# Load data ---------------------------------------------------------------

umap_embeddings    <- readRDS(here("data", "processed", "umap_embeddings.rds"))
metadata_processed <- readRDS(here("data", "processed", "metadata_processed.rds"))
elbow_df           <- readRDS(here("data", "processed", "elbow_df.rds"))
min_pc             <- unique(elbow_df$min_pc)


# Wrangle data ------------------------------------------------------------

umap_df <- data.frame(
  UMAP_1     = umap_embeddings[, 1],
  UMAP_2     = umap_embeddings[, 2],
  sample_id  = metadata_processed$orig.ident,
  sex        = metadata_processed$sex,
  cluster    = as.character(metadata_processed[[cluster_col]]),
  percent_mt = metadata_processed$percent.mt
)

umap_df$cell_type <- cluster_annotations[umap_df$cluster]
umap_df$cell_type <- factor(umap_df$cell_type, levels = names(celltype_colors))

cluster_centroids <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ cluster,
  data = umap_df,
  FUN  = median
)

celltype_centroids <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ cell_type,
  data = umap_df,
  FUN  = median
)

## Arrow start positions from data range
umap1_min <- min(umap_df$UMAP_1)
umap2_min <- min(umap_df$UMAP_2)


# Helper functions --------------------------------------------------------

## Minimal arrow axes — Scavuzzo et al. style
umap_arrows <- function(x_start  = umap1_min,
                        y_start  = umap2_min,
                        length   = arrow_length,
                        txt_size = 2.5) {
  list(
    annotate(
      geom = "segment",
      x = x_start, xend = x_start + length,
      y = y_start, yend = y_start,
      arrow     = arrow(length = unit(0.15, "cm"), type = "closed"),
      color     = "black",
      linewidth = 0.4
    ),
    annotate(
      geom = "segment",
      x = x_start, xend = x_start,
      y = y_start, yend = y_start + length,
      arrow     = arrow(length = unit(0.15, "cm"), type = "closed"),
      color     = "black",
      linewidth = 0.4
    ),
    annotate(
      geom  = "text",
      x     = x_start + length / 2,
      y     = y_start - 0.8,
      label = "UMAP 1",
      size  = txt_size,
      hjust = 0.5,
      color = "black"
    ),
    annotate(
      geom  = "text",
      x     = x_start - 0.8,
      y     = y_start + length / 2,
      label = "UMAP 2",
      size  = txt_size,
      angle = 90,
      hjust = 0.5,
      color = "black"
    )
  )
}


# Visualization -----------------------------------------------------------

## Shared minimal theme — no titles, captions, axes, or border
umap_theme <- theme_bw(base_size = 9) +
  theme(
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    plot.title   = element_blank(),
    plot.caption = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    axis.line    = element_blank(),
    legend.key.size = unit(0.4, "cm")
  )

## Elbow plot — keep standard axes for readability
p_elbow <- elbow_df |>
  ggplot(aes(x = cumu, y = pct, label = rank, color = rank > min_pc)) +
  geom_text(size = 2.5) +
  geom_vline(xintercept = 90, color = "grey60", linetype = "dashed") +
  geom_hline(
    yintercept = min(elbow_df$pct[elbow_df$pct > 5]),
    color      = "grey60",
    linetype   = "dashed"
  ) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "grey70"),
    labels = c("FALSE" = "Selected PCs", "TRUE" = "Excluded PCs")
  ) +
  labs(x = "Cumulative % Variance", y = "% Variance per PC", color = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_blank())

## UMAP colored by sample
p_umap_sample <- umap_df |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = sample_id)) +
  geom_point(size = 0.2, alpha = 0.5) +
  umap_arrows() +
  labs(color = "Sample") +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3)))

## UMAP colored by sex
p_umap_sex <- umap_df |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = sex)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("Female" = "hotpink3", "Male" = "steelblue")) +
  umap_arrows() +
  labs(color = "Sex") +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3)))

## UMAP colored by cluster with centroid labels
p_umap_clusters <- umap_df |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(
    data        = cluster_centroids,
    aes(x = UMAP_1, y = UMAP_2, label = cluster),
    color       = "black",
    size        = 2.5,
    fontface    = "bold",
    inherit.aes = FALSE
  ) +
  umap_arrows() +
  labs(color = "Cluster") +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))

## UMAP colored by percent mitochondrial reads
p_umap_pct_mt <- umap_df |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = percent_mt)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_color_viridis_c(option = "magma", name = "% MT") +
  umap_arrows() +
  umap_theme

## Annotated UMAP colored by cell type
p_umap_annotated <- umap_df |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_color_manual(values = celltype_colors, name = "Cell Type") +
  umap_arrows() +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3)))

## Annotated UMAP with cell type labels at centroids
p_umap_annotated_labeled <- umap_df |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.4) +
  scale_color_manual(values = celltype_colors, name = "Cell Type") +
  geom_label(
    data        = celltype_centroids,
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    color       = "black",
    fill        = "white",
    size        = 2.5,
    fontface    = "bold",
    label.size  = 0.2,
    alpha       = 0.8,
    inherit.aes = FALSE
  ) +
  umap_arrows() +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Save outputs ------------------------------------------------------------

ggsave(here("plots", "pca_elbow_plot.pdf"),        p_elbow,                 width = 7,  height = 5)
ggsave(here("plots", "umap_by_sample.pdf"),         p_umap_sample,           width = 8,  height = 6)
ggsave(here("plots", "umap_by_sex.pdf"),            p_umap_sex,              width = 7,  height = 6)
ggsave(here("plots", "umap_clusters.pdf"),          p_umap_clusters,         width = 9,  height = 7)
ggsave(here("plots", "umap_pct_mt.pdf"),            p_umap_pct_mt,           width = 7,  height = 6)
ggsave(here("plots", "umap_annotated.pdf"),         p_umap_annotated,        width = 9,  height = 7)
ggsave(here("plots", "umap_annotated_labeled.pdf"), p_umap_annotated_labeled, width = 10, height = 7)

message("UMAP plots saved to plots/")


# Session info ------------------------------------------------------------

sessionInfo()