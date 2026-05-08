# Marker Genes ------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Finds marker genes for each cluster in the processed Seurat
#              object using FindAllMarkers, and summarizes the top 10 markers
#              per cluster by log2 fold change for downstream visualization.
#              PrepSCTFindMarkers() is required before FindAllMarkers() when
#              SCTransform was run per-sample, to reconcile unequal library
#              sizes across SCT models (Seurat v5).
#              Set input_suffix to "_denoised" for the cellsweep pass.
# Input:       data/processed/seurat_processed<input_suffix>.rds
# Output:      data/processed/seurat_markers<input_suffix>.rds
#              data/processed/top10_markers<input_suffix>.rds
#              data/processed/cluster_markers<input_suffix>.txt
# Note:        Submit as a SLURM job — PrepSCTFindMarkers() and
#              FindAllMarkers() are both memory and time intensive.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Input/output suffix — set to "_denoised" for cellsweep pass
input_suffix <- "_denoised"

## Default clustering resolution — must match normalize_cluster.R
default_resolution <- 0.5

## Minimum fraction of cells expressing a gene to test
min_pct <- 0.25

## Minimum log2 fold change threshold for marker testing
logfc_threshold <- 0.25


# Libraries ---------------------------------------------------------------

library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

seurat_path <- here("data", "processed", paste0("seurat_processed", input_suffix, ".rds"))

message("Loading ", basename(seurat_path), "...")
seurat_processed <- readRDS(seurat_path)
message("Loaded: ", ncol(seurat_processed), " cells, ", nrow(seurat_processed), " features")


# Analysis ----------------------------------------------------------------

## Set cluster identity explicitly before finding markers
Idents(seurat_processed) <- paste0("SCT_snn_res.", default_resolution)

message("Using resolution ", default_resolution, " — ",
        length(levels(Idents(seurat_processed))), " clusters")

## Required when SCTransform was run per-sample (Seurat v5)
message("Running PrepSCTFindMarkers — this may take a while...")
seurat_processed <- PrepSCTFindMarkers(seurat_processed)
message("PrepSCTFindMarkers complete")

## Find marker genes for each cluster vs all other cells
message("Running FindAllMarkers...")
seurat_markers <- FindAllMarkers(
  seurat_processed,
  only.pos        = TRUE,
  min.pct         = min_pct,
  logfc.threshold = logfc_threshold,
  verbose         = FALSE
)

message("Found ", nrow(seurat_markers), " marker genes across ",
        length(unique(seurat_markers$cluster)), " clusters")

## Top 10 markers per cluster ordered by log2 fold change
top10_markers <- seurat_markers |>
  (\(df) split(df, df$cluster))() |>
  lapply(\(x) x[order(-x$avg_log2FC), ][1:min(10, nrow(x)), ]) |>
  do.call(what = rbind)


# Save outputs ------------------------------------------------------------

saveRDS(
  seurat_markers,
  file = here("data", "processed", paste0("seurat_markers", input_suffix, ".rds"))
)
saveRDS(
  top10_markers,
  file = here("data", "processed", paste0("top10_markers", input_suffix, ".rds"))
)
write.table(
  seurat_markers,
  file      = here("data", "processed", paste0("cluster_markers", input_suffix, ".txt")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("Saved seurat_markers", input_suffix, ".rds, top10_markers", input_suffix,
        ".rds, cluster_markers", input_suffix, ".txt")


# Session info ------------------------------------------------------------

sessionInfo()