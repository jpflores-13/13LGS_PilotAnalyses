# Extract Lightweight Data ------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: One-time extraction script that loads the full 20 GB Seurat
#              object and saves lightweight files needed by the plot scripts.
#              Extracts UMAP embeddings, per-cell metadata, and SCT expression
#              for ALL significant marker genes plus known ileum marker genes.
#              Run this as a SLURM job once, then iterate on visualizations
#              interactively without ever loading the full object again.
# Input:       data/processed/seurat_processed.rds
#              data/processed/seurat_markers.rds
# Output:      data/processed/umap_embeddings.rds    — UMAP coordinates
#              data/processed/metadata_processed.rds  — per-cell metadata
#              data/processed/marker_expr.rds         — expression matrix for
#                                                       all marker genes +
#                                                       known ileum markers
# Note:        Submit as a SLURM job — loads the full 20 GB Seurat object.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

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

## Load full Seurat object — memory intensive, run as SLURM job
message("Loading seurat_processed.rds...")
seurat_processed <- readRDS(here("data", "processed", "seurat_processed.rds"))
message("Loaded: ", ncol(seurat_processed), " cells, ", nrow(seurat_processed), " features")

## Load marker gene table to get all significant marker genes
seurat_markers <- readRDS(here("data", "processed", "seurat_markers.rds"))
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
  file = here("data", "processed", "umap_embeddings.rds")
)
saveRDS(
  metadata_processed,
  file = here("data", "processed", "metadata_processed.rds")
)
saveRDS(
  marker_expr,
  file = here("data", "processed", "marker_expr.rds")
)

message("Lightweight extracts saved to data/processed/")
message("  umap_embeddings.rds  — ", nrow(umap_embeddings), " cells x 2 UMAP dimensions")
message("  metadata_processed.rds — ", nrow(metadata_processed), " cells")
message("  marker_expr.rds      — ", nrow(marker_expr), " cells x ", ncol(marker_expr), " genes")
message("Plot scripts can now be run interactively without loading the full object")


# Session info ------------------------------------------------------------

sessionInfo()