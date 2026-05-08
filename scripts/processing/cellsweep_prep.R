# CellSweep Prep -----------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-08
# Project:     13LGS_PilotAnalyses
# Description: Prepares per-sample input for cellsweep denoising. Reads raw
#              counts from filtered and unfiltered CellRanger .h5 files,
#              maps first-pass cluster annotations from the processed Seurat
#              object onto filtered barcodes, and writes per-sample flat
#              files (MTX + metadata) for conversion to .h5ad by
#              cellsweep_convert.py. No reticulate dependency.
#
#              The output mirrors the cellsweep tutorial structure:
#                - Filtered (real) cells get a celltype label from cluster
#                  annotations; is_empty = FALSE
#                - Empty droplets (barcodes in raw but not filtered) get
#                  celltype = "empty"; is_empty = TRUE
#
# Input:       data/raw/<sample_id>_filtered_feature_bc_matrix.h5
#              data/raw/<sample_id>_raw_feature_bc_matrix.h5
#              data/processed/metadata_processed.rds
# Output:      data/processed/cellsweep/<sample_id>/
#                counts.mtx          — sparse count matrix (cells+empty x genes)
#                barcodes.tsv        — all barcodes (filtered + empty)
#                genes.tsv           — gene names
#                metadata.csv        — per-barcode metadata with celltype + is_empty
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Default clustering resolution — must match normalize_cluster.R
default_resolution <- 0.5
cluster_col        <- paste0("SCT_snn_res.", default_resolution)

## Directory containing raw .h5 files
h5_dir <- here::here("data", "raw")

## Suffix patterns for h5 files
filtered_suffix   <- "_filtered_feature_bc_matrix.h5"
unfiltered_suffix <- "_raw_feature_bc_matrix.h5"

## First-pass cluster annotations — must match umap_plots.R
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


# Libraries ---------------------------------------------------------------

library(here)
library(Matrix)
library(Seurat)


# Load data ---------------------------------------------------------------

metadata_processed <- readRDS(here("data", "processed", "metadata_processed.rds"))
sample_ids         <- unique(metadata_processed$orig.ident)

message("Preparing cellsweep input for ", length(sample_ids), " samples")

dir.create(here("data", "processed", "cellsweep"), showWarnings = FALSE)


# Process each sample -----------------------------------------------------

for (sample_id in sample_ids) {
  
  message("\nProcessing: ", sample_id)
  
  out_dir <- here("data", "processed", "cellsweep", sample_id)
  dir.create(out_dir, showWarnings = FALSE)
  
  ## --- Filtered cells ---------------------------------------------------
  
  filtered_path <- file.path(h5_dir, paste0(sample_id, filtered_suffix))
  
  if (!file.exists(filtered_path)) {
    message("  WARNING: filtered h5 not found — skipping: ", filtered_path)
    next
  }
  
  message("  Reading filtered h5...")
  counts_filtered <- Read10X_h5(filtered_path)
  
  ## Get metadata for this sample
  sample_meta     <- metadata_processed[metadata_processed$orig.ident == sample_id, ]
  seurat_barcodes <- rownames(sample_meta)
  
  ## Strip all prefixes — extract the raw barcode (ACGT...-1) regardless of
  ## how many sample prefixes were stacked during Seurat merge
  stripped_barcodes <- gsub("^.*_([ACGT]+-[0-9]+)$", "\\1", seurat_barcodes)
  
  ## Match stripped Seurat barcodes to h5 barcodes
  matched        <- stripped_barcodes %in% colnames(counts_filtered)
  stripped_match <- stripped_barcodes[matched]
  sample_meta    <- sample_meta[matched, ]
  
  message("  Cells in Seurat:     ", length(seurat_barcodes))
  message("  Cells matched to h5: ", sum(matched))
  
  if (sum(matched) == 0) {
    message("  WARNING: No barcodes matched — skipping")
    next
  }
  
  ## Subset filtered counts to matched cells
  counts_matched <- counts_filtered[, stripped_match, drop = FALSE]
  
  ## Map cluster IDs to cell type labels
  cluster_ids <- as.character(sample_meta[[cluster_col]])
  cell_types  <- cluster_annotations[cluster_ids]
  
  ## Build filtered cell metadata
  meta_filtered <- data.frame(
    barcode  = stripped_match,
    celltype = cell_types,
    is_empty = FALSE,
    row.names = stripped_match,
    stringsAsFactors = FALSE
  )
  
  ## --- Empty droplets ---------------------------------------------------
  
  unfiltered_path <- file.path(h5_dir, paste0(sample_id, unfiltered_suffix))
  
  if (!file.exists(unfiltered_path)) {
    message("  WARNING: unfiltered h5 not found — cannot include empty droplets")
    counts_all   <- counts_matched
    meta_all     <- meta_filtered
    all_barcodes <- stripped_match
  } else {
    message("  Reading unfiltered h5...")
    counts_unfiltered <- Read10X_h5(unfiltered_path)
    
    ## Empty droplets = barcodes in raw but NOT in filtered h5
    ## Use full filtered h5 barcodes (not just Seurat-matched) to avoid
    ## misclassifying QC-filtered real cells as empty droplets
    empty_barcodes <- setdiff(colnames(counts_unfiltered), colnames(counts_filtered))
    
    message("  Empty droplets: ", length(empty_barcodes))
    
    counts_empty <- counts_unfiltered[, empty_barcodes, drop = FALSE]
    
    meta_empty <- data.frame(
      barcode  = empty_barcodes,
      celltype = "empty",
      is_empty = TRUE,
      row.names = empty_barcodes,
      stringsAsFactors = FALSE
    )
    
    ## Combine filtered cells + empty droplets
    shared_genes  <- intersect(rownames(counts_matched), rownames(counts_empty))
    counts_all    <- cbind(
      counts_matched[shared_genes, , drop = FALSE],
      counts_empty[shared_genes, , drop = FALSE]
    )
    meta_all     <- rbind(meta_filtered, meta_empty)
    all_barcodes <- c(stripped_match, empty_barcodes)
  }
  
  ## --- Write outputs ----------------------------------------------------
  
  ## cellsweep expects cells x genes — transpose before writing
  counts_t <- t(counts_all)
  
  Matrix::writeMM(counts_t, file.path(out_dir, "counts.mtx"))
  
  write.table(
    all_barcodes,
    file      = file.path(out_dir, "barcodes.tsv"),
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
  write.table(
    rownames(counts_all),
    file      = file.path(out_dir, "genes.tsv"),
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
  write.csv(meta_all, file.path(out_dir, "metadata.csv"), row.names = FALSE)
  
  message("  Cells written:          ", sum(!meta_all$is_empty))
  message("  Empty droplets written: ", sum(meta_all$is_empty))
  message("  Done: ", out_dir)
}

message("\nCellSweep prep complete!")
message("Next step: run cellsweep_convert.py to build per-sample .h5ad files")


# Session info ------------------------------------------------------------

sessionInfo()