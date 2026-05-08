# Extract Lightweight Data ------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: One-time extraction script that loads the full Seurat object
#              and saves lightweight files needed by the plot scripts.
#              Extracts UMAP embeddings, per-cell metadata, and SCT expression
#              for ALL significant marker genes plus known ileum marker genes.
#              Run this as a SLURM job once, then iterate on visualizations
#              interactively without ever loading the full object again.
#              Set input_suffix to "_denoised" to extract from the cellsweep
#              denoised object.
# Input:       data/processed/seurat_processed<input_suffix>.rds
#              data/processed/seurat_markers<input_suffix>.rds
# Output:      data/processed/umap_embeddings<input_suffix>.rds
#              data/processed/metadata_processed<input_suffix>.rds
#              data/processed/marker_expr<input_suffix>.rds
# Note:        Submit as a SLURM job — loads the full Seurat object.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Input/output suffix — set to "_denoised" for cellsweep pass
input_suffix <- "_denoised"

## Known ileum marker genes for dot plot
## Must match gene_groups in marker_plots.R
ileum_markers <- c(
  "FABP1", "FABP2", "APOA1", "SLC5A1",
  "MUC2", "TFF3", "CLCA1",
  "DEFA5", "DEFA6", "LYZ",
  "CHGA", "CHGB", "SCG2",
  "DCLK1", "POU2F3", "TRPM5",
  "LGR5", "OLFM4", "ASCL2",
  "SOX10", "GFAP", "PLP1", "S100B",
  "ELAVL4", "TUBB3",
  "ACTA2", "MYH11",
  "PECAM1", "CDH5",
  "CD68", "CSF1R", "PTPRC",
  "VIM", "COL1A1", "PDGFRA"
)


# Libraries ---------------------------------------------------------------

library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

seurat_path  <- here("data", "processed", paste0("seurat_processed", input_suffix, ".rds"))
markers_path <- here("data", "processed", paste0("seurat_markers", input_suffix, ".rds"))

message("Loading ", basename(seurat_path), "...")
seurat_processed <- readRDS(seurat_path)
message("Loaded: ", ncol(seurat_processed), " cells, ", nrow(seurat_processed), " features")

seurat_markers <- readRDS(markers_path)
message("Loaded seurat_markers: ", nrow(seurat_markers), " marker genes")


# Wrangle data ------------------------------------------------------------

## Extract UMAP embeddings
umap_embeddings <- seurat_processed@reductions$umap@cell.embeddings
message("Extracted UMAP embeddings: ", nrow(umap_embeddings), " cells")

## Extract per-cell metadata
metadata_processed <- seurat_processed@meta.data
message("Extracted metadata: ", nrow(metadata_processed), " cells")

## Build combined gene list:
## ALL significant marker genes + known ileum markers for dot plot
all_marker_genes <- unique(seurat_markers$gene)
genes_to_extract <- unique(c(ileum_markers, all_marker_genes))

## Filter to genes actually present in the feature matrix
genes_to_extract <- genes_to_extract[genes_to_extract %in% rownames(seurat_processed)]

message(
  "Extracting SCT expression for ", length(genes_to_extract), " genes",
  " (", length(all_marker_genes), " marker genes + ",
  sum(ileum_markers %in% rownames(seurat_processed)), " ileum markers detected)"
)

## Extract SCT normalized expression — rows = cells, cols = genes
marker_expr <- t(
  GetAssayData(seurat_processed, assay = "SCT", layer = "data")[genes_to_extract, ]
)

message("Expression matrix: ", nrow(marker_expr), " cells x ", ncol(marker_expr), " genes")


# Save outputs ------------------------------------------------------------

saveRDS(
  umap_embeddings,
  file = here("data", "processed", paste0("umap_embeddings", input_suffix, ".rds"))
)
saveRDS(
  metadata_processed,
  file = here("data", "processed", paste0("metadata_processed", input_suffix, ".rds"))
)
saveRDS(
  marker_expr,
  file = here("data", "processed", paste0("marker_expr", input_suffix, ".rds"))
)

message("Lightweight extracts saved to data/processed/")
message("  umap_embeddings", input_suffix, ".rds  — ", nrow(umap_embeddings), " cells x 2 UMAP dimensions")
message("  metadata_processed", input_suffix, ".rds — ", nrow(metadata_processed), " cells")
message("  marker_expr", input_suffix, ".rds      — ", nrow(marker_expr), " cells x ", ncol(marker_expr), " genes")
message("Plot scripts can now be run interactively without loading the full object")


# Session info ------------------------------------------------------------

sessionInfo()