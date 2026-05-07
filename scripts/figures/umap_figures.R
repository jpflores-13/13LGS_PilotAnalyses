# UMAP Figures ------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Generates UMAP visualizations from the processed and
#              cellsweep-denoised Seurat objects. Produces UMAPs colored
#              by sample, sex, and cluster identity, plus an elbow plot
#              showing the selected number of PCs.
# Input:       data/processed/seurat_processed.rds
#              data/processed/seurat_cellsweep.rds
#              data/processed/elbow_df.rds
# Output:      plots/pca_elbow_plot.pdf
#              plots/umap_by_sample.pdf
#              plots/umap_by_sex.pdf
#              plots/umap_clusters.pdf
#              plots/umap_cellsweep.pdf
# Note:        This project uses renv for reproducibility.
#              Run renv::restore() before executing this script.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Default clustering resolution — must match what was set in normalize_cluster.R
default_resolution <- 0.5


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Processed Seurat object from normalize_cluster.R
seurat_processed <- readRDS(here("data", "processed", "seurat_processed.rds"))

## Cellsweep-denoised Seurat object from cellsweep_benchmark.R
seurat_cellsweep <- readRDS(here("data", "processed", "seurat_cellsweep.rds"))

## PC selection data frame from normalize_cluster.R
elbow_df <- readRDS(here("data", "processed", "elbow_df.rds"))
min_pc   <- unique(elbow_df$min_pc)


# Visualization -----------------------------------------------------------

## Elbow plot with selected PC cutoff highlighted
p_elbow <- elbow_df |>
  ggplot(aes(x = cumu, y = pct, label = rank, color = rank > min_pc)) +
  geom_text() +
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
  labs(
    x     = "Cumulative % Variance",
    y     = "% Variance per PC",
    title = paste0("PC Selection — using first ", min_pc, " PCs"),
    color = NULL
  ) +
  theme_bw()

## UMAP colored by sample identity
p_umap_sample <- DimPlot(
  seurat_processed,
  reduction = "umap",
  group.by  = "orig.ident",
  pt.size   = 0.5
) + labs(title = "UMAP — colored by sample")

## UMAP colored by sex
p_umap_sex <- DimPlot(
  seurat_processed,
  reduction = "umap",
  group.by  = "sex",
  cols      = c("Female" = "hotpink3", "Male" = "steelblue"),
  pt.size   = 0.5
) + labs(title = "UMAP — colored by sex")

## UMAP colored by cluster identity
p_umap_clusters <- DimPlot(
  seurat_processed,
  reduction = "umap",
  pt.size   = 0.5
) + labs(title = paste0("UMAP — clusters (res = ", default_resolution, ")"))

## UMAP of cellsweep-denoised object for benchmarking comparison
p_umap_cellsweep <- DimPlot(
  seurat_cellsweep,
  reduction = "umap",
  pt.size   = 0.5
) + labs(title = paste0("UMAP — cellsweep denoised (res = ", default_resolution, ")"))


# Save outputs ------------------------------------------------------------

ggsave(
  filename = here("plots", "pca_elbow_plot.pdf"),
  plot     = p_elbow,
  width    = 7,
  height   = 5
)
ggsave(
  filename = here("plots", "umap_by_sample.pdf"),
  plot     = p_umap_sample,
  width    = 8,
  height   = 6
)
ggsave(
  filename = here("plots", "umap_by_sex.pdf"),
  plot     = p_umap_sex,
  width    = 7,
  height   = 6
)
ggsave(
  filename = here("plots", "umap_clusters.pdf"),
  plot     = p_umap_clusters,
  width    = 7,
  height   = 6
)
ggsave(
  filename = here("plots", "umap_cellsweep.pdf"),
  plot     = p_umap_cellsweep,
  width    = 7,
  height   = 6
)


# Session info ------------------------------------------------------------

sessionInfo()