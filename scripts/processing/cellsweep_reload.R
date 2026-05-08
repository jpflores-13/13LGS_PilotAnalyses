# CellSweep Reload --------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-08
# Project:     13LGS_PilotAnalyses
# Description: Loads per-sample denoised count matrices from cellsweep
#              (_denoised.h5ad) and injects them as the RNA count layer in
#              the corresponding doublet-cleaned Seurat objects. Produces
#              one Seurat object per sample with clean counts, ready for
#              re-running sctransform_per_sample.R.
#
#              Uses zellkonverter to read .h5ad — no reticulate dependency.
#
# Input:       data/processed/cellsweep/<sample_id>_denoised.h5ad
#              data/processed/doublets/<sample_id>_clean.rds
# Output:      data/processed/cellsweep_clean/<sample_id>_clean.rds
# Note:        Run after all cellsweep_run.py jobs complete.
#              Then resubmit sctransform_per_sample.R pointed at
#              data/processed/cellsweep_clean/ instead of doublets/.
# -------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------

library(here)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)


# Setup -------------------------------------------------------------------

denoised_dir <- here("data", "processed", "cellsweep")
doublet_dir  <- here("data", "processed", "doublets")
output_dir   <- here("data", "processed", "cellsweep_clean")

dir.create(output_dir, showWarnings = FALSE)

denoised_files <- list.files(
  path       = denoised_dir,
  pattern    = "_denoised\\.h5ad$",
  full.names = TRUE
)

if (length(denoised_files) == 0) {
  stop("No _denoised.h5ad files found in ", denoised_dir)
}

sample_ids <- gsub("_denoised\\.h5ad$", "", basename(denoised_files))

message("Found ", length(denoised_files), " denoised samples: ",
        paste(sample_ids, collapse = ", "))


# Process each sample -----------------------------------------------------

for (i in seq_along(sample_ids)) {
  
  sample_id <- sample_ids[i]
  message("\nProcessing: ", sample_id)
  
  ## Load denoised AnnData as SingleCellExperiment
  ## zellkonverter reads .X as counts assay, obs as colData
  sce <- zellkonverter::readH5AD(denoised_files[i])
  
  ## Extract denoised count matrix — genes x cells (Seurat convention)
  ## zellkonverter transposes AnnData's cells x genes to genes x cells
  counts_denoised <- assay(sce, "X")
  
  message("  Denoised matrix: ", nrow(counts_denoised), " genes x ",
          ncol(counts_denoised), " cells")
  
  ## Load matching doublet-cleaned Seurat object
  clean_path <- file.path(doublet_dir, paste0(sample_id, "_clean.rds"))
  
  if (!file.exists(clean_path)) {
    message("  WARNING: ", clean_path, " not found — skipping")
    next
  }
  
  seurat_obj <- readRDS(clean_path)
  
  ## Strip Seurat sample prefix to match raw h5 barcodes used by cellsweep
  seurat_barcodes   <- colnames(seurat_obj)
  stripped_barcodes <- gsub(paste0("^", sample_id, "_"), "", seurat_barcodes)
  
  ## Find common barcodes between denoised output and Seurat object
  denoised_barcodes <- colnames(counts_denoised)
  common_barcodes   <- intersect(stripped_barcodes, denoised_barcodes)
  
  message("  Seurat cells:   ", ncol(seurat_obj))
  message("  Denoised cells: ", ncol(counts_denoised))
  message("  Matched cells:  ", length(common_barcodes))
  
  if (length(common_barcodes) == 0) {
    message("  WARNING: No barcodes matched — skipping")
    next
  }
  
  ## Subset Seurat to matched barcodes (preserves metadata, reductions, etc.)
  seurat_idx <- match(common_barcodes, stripped_barcodes)
  seurat_sub <- seurat_obj[, seurat_idx]
  
  ## Align genes between denoised matrix and Seurat feature set
  common_genes    <- intersect(rownames(counts_denoised), rownames(seurat_sub))
  counts_denoised <- counts_denoised[common_genes, common_barcodes, drop = FALSE]
  seurat_sub      <- seurat_sub[common_genes, ]
  
  message("  Common genes: ", length(common_genes))
  
  ## Replace RNA count layer with denoised counts
  seurat_sub <- SetAssayData(
    object   = seurat_sub,
    assay    = "RNA",
    layer    = "counts",
    new.data = counts_denoised
  )
  
  ## Save
  out_path <- file.path(output_dir, paste0(sample_id, "_clean.rds"))
  saveRDS(seurat_sub, file = out_path)
  message("  Saved: ", out_path)
}

message("\nCellSweep reload complete!")
message("Next step: resubmit sctransform_per_sample.R reading from data/processed/cellsweep_clean/")


# Session info ------------------------------------------------------------

sessionInfo()