# Marker Genes ------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Finds marker genes for each cluster in the processed Seurat
#              object using FindAllMarkers, and summarizes the top 10 markers
#              per cluster by log2 fold change for downstream visualization.
# Input:       data/processed/seurat_processed.rds
# Output:      data/processed/seurat_markers.rds   — full marker gene table
#              data/processed/top10_markers.rds     — top 10 per cluster
#              data/processed/cluster_markers.txt   — tab-separated marker table
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Minimum fraction of cells expressing a gene to test
min_pct <- 0.25

## Minimum log2 fold change threshold for marker testing
logfc_threshold <- 0.25


# Libraries ---------------------------------------------------------------

library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Processed Seurat object from normalize_cluster.R
seurat_processed <- readRDS(here("data", "processed", "seurat_processed.rds"))


# Analysis ----------------------------------------------------------------

## Find marker genes for each cluster vs all other cells
## only.pos = TRUE: report only upregulated markers per cluster
seurat_markers <- FindAllMarkers(
  seurat_processed,
  only.pos        = TRUE,
  min.pct         = min_pct,
  logfc.threshold = logfc_threshold,
  verbose         = FALSE
)

## Top 10 markers per cluster ordered by log2 fold change
## Used downstream for heatmap visualization
top10_markers <- seurat_markers |>
  (\(df) split(df, df$cluster))() |>
  lapply(\(x) x[order(-x$avg_log2FC), ][1:min(10, nrow(x)), ]) |>
  do.call(what = rbind)


# Save outputs ------------------------------------------------------------

saveRDS(seurat_markers, file = here("data", "processed", "seurat_markers.rds"))
saveRDS(top10_markers,  file = here("data", "processed", "top10_markers.rds"))

write.table(
  seurat_markers,
  file      = here("data", "processed", "cluster_markers.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)


# Session info ------------------------------------------------------------

sessionInfo()