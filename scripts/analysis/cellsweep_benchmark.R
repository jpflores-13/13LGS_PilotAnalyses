# Cellsweep Benchmarking --------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Benchmarks cellsweep denoising against the standard Seurat
#              pipeline. Converts the processed Seurat object to AnnData,
#              runs cellsweep using cluster labels as proxy cell type
#              annotations, pulls denoised counts back into R, and runs
#              an identical Seurat pipeline for like-for-like UMAP comparison.
# Input:       data/processed/seurat_processed.rds
# Output:      data/processed/seurat_cellsweep.rds  — denoised Seurat object
#              data/processed/seurat_merged_raw.h5ad — AnnData input for cellsweep
#              data/processed/seurat_cellsweep.h5ad  — cellsweep denoised output
# Note:        Requires Python with cellsweep installed: pip install cellsweep
#              Configure reticulate to point at the correct Python environment
#              before running: reticulate::use_virtualenv() or use_condaenv()
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Variables to regress out during SCTransform (must match normalize_cluster.R)
vars_to_regress <- "percent.mt"

## Harmony batch correction grouping variables (must match normalize_cluster.R)
harmony_vars <- c("orig.ident", "sex")


# Libraries ---------------------------------------------------------------

library(here)
library(reticulate)
library(Seurat)


# Load data ---------------------------------------------------------------

## Processed Seurat object from normalize_cluster.R
## Cluster labels from this object are used as proxy cell types for cellsweep
seurat_processed <- readRDS(here("data", "processed", "seurat_processed.rds"))

## Recover min_pc from elbow_df saved in normalize_cluster.R
elbow_df <- readRDS(here("data", "processed", "elbow_df.rds"))
min_pc   <- unique(elbow_df$min_pc)


# Analysis ----------------------------------------------------------------

## Import Python modules via reticulate
cellsweep_py <- import("cellsweep")
ad           <- import("anndata")

## Extract raw counts from the RNA assay for cellsweep input
## Cellsweep expects raw (unnormalized) counts, not SCTransform-normalized
counts_mat <- GetAssayData(seurat_processed, assay = "RNA", layer = "counts")

## Build AnnData object from raw counts and Seurat metadata
adata_raw <- ad$AnnData(
  X   = Matrix::t(counts_mat),
  obs = seurat_processed@meta.data
)

## Add cluster labels as proxy cell type annotations
## Replace with biological cell type labels once annotation is complete
adata_raw$obs[["celltype"]] <- as.character(Idents(seurat_processed))

## Write AnnData to disk — cellsweep reads from .h5ad
h5ad_raw_path      <- here("data", "processed", "seurat_merged_raw.h5ad")
h5ad_denoised_path <- here("data", "processed", "seurat_cellsweep.h5ad")
adata_raw$write_h5ad(h5ad_raw_path)

## Run cellsweep denoising
message("Running cellsweep denoising...")
adata_cellsweep <- cellsweep_py$denoise_count_matrix(
  h5ad_raw_path,
  adata_out = h5ad_denoised_path
)

## Pull denoised counts back into R and restore row/col names
counts_denoised <- Matrix::t(
  as(adata_cellsweep$X, "CsparseMatrix")
)
rownames(counts_denoised) <- rownames(counts_mat)
colnames(counts_denoised) <- colnames(counts_mat)

## Build a parallel Seurat object from denoised counts
seurat_cellsweep <- CreateSeuratObject(
  counts    = counts_denoised,
  meta.data = seurat_processed@meta.data
)

## Run identical pipeline on denoised counts for like-for-like comparison
message("Running Seurat pipeline on cellsweep-denoised counts...")

seurat_cellsweep <- SCTransform(
  seurat_cellsweep,
  vars.to.regress       = vars_to_regress,
  return.only.var.genes = FALSE,
  verbose               = FALSE
)
seurat_cellsweep <- RunPCA(seurat_cellsweep, verbose = FALSE)
seurat_cellsweep <- RunHarmony(
  seurat_cellsweep,
  group.by.vars = harmony_vars,
  assay.use     = "SCT",
  verbose       = FALSE
)
seurat_cellsweep <- RunUMAP(
  seurat_cellsweep,
  reduction = "harmony",
  dims      = 1:min_pc,
  verbose   = FALSE
)
seurat_cellsweep <- FindNeighbors(
  seurat_cellsweep,
  reduction = "harmony",
  dims      = 1:min_pc,
  verbose   = FALSE
)
seurat_cellsweep <- FindClusters(
  seurat_cellsweep,
  resolution = 0.5,
  verbose    = FALSE
)


# Save outputs ------------------------------------------------------------

saveRDS(seurat_cellsweep, file = here("data", "processed", "seurat_cellsweep.rds"))


# Session info ------------------------------------------------------------

sessionInfo()