# CellSweep Prep ----------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-08
# Project:     13LGS_PilotAnalyses
# Description: Prepares per-sample data for cellsweep denoising. Reads raw
#              counts from filtered and unfiltered .h5 files, adds cell type
#              labels from first-pass annotation, and exports sparse count
#              matrices and metadata as flat files for conversion to .h5ad
#              by cellsweep_convert.py. No reticulate dependency.
# Input:       data/raw/<sample_id>_filtered_feature_bc_matrix.h5
#              data/raw/<sample_id>_raw_feature_bc_matrix.h5 (if available)
#              data/processed/metadata_processed.rds
# Output:      data/processed/cellsweep/<sample_id>/
#                  counts_filtered.mtx    — sparse count matrix (cells x genes)
#                  barcodes_filtered.tsv  — cell barcodes
#                  genes.tsv              — gene names
#                  metadata.csv           — per-cell metadata + cell type labels
#                  counts_empty.mtx       — empty droplet counts (if unfiltered)
#                  barcodes_empty.tsv     — empty droplet barcodes (if unfiltered)
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Default clustering resolution — must match normalize_cluster.R
default_resolution <- 0.5
cluster_col        <- paste0("SCT_snn_res.", default_resolution)

## Directory containing raw .h5 files
filtered_h5_dir   <- here::here("data", "raw")
unfiltered_h5_dir <- here::here("data", "raw")

## Set to TRUE if unfiltered .h5 files are available
has_unfiltered <- TRUE

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
  
  ## Create per-sample output directory
  out_dir <- here("data", "processed", "cellsweep", sample_id)
  dir.create(out_dir, showWarnings = FALSE)
  
  ## Path to filtered h5
  h5_path <- file.path(filtered_h5_dir, paste0(sample_id, filtered_suffix))
  
  if (!file.exists(h5_path)) {
    message("  WARNING: h5 not found — skipping: ", h5_path)
    next
  }
  
  ## Read raw counts from filtered h5
  message("  Reading filtered h5...")
  counts_filtered <- Read10X_h5(h5_path)
  
  ## Get metadata for this sample
  sample_meta     <- metadata_processed[metadata_processed$orig.ident == sample_id, ]
  seurat_barcodes <- rownames(sample_meta)
  
  ## Strip sample prefix from Seurat barcodes to match raw h5 barcodes
  stripped_barcodes <- gsub(paste0("^", sample_id, "_"), "", seurat_barcodes)
  matched           <- stripped_barcodes %in% colnames(counts_filtered)
  
  message("  Cells in metadata:   ", nrow(sample_meta))
  message("  Cells matched to h5: ", sum(matched))
  
  if (sum(matched) == 0) {
    message("  WARNING: No barcodes matched — skipping")
    next
  }
  
  ## Subset to matched cells
  sample_meta_matched <- sample_meta[matched, ]
  stripped_matched    <- stripped_barcodes[matched]
  counts_matched      <- counts_filtered[, stripped_matched]
  
  ## Add cell type labels
  cluster_ids <- as.character(sample_meta_matched[[cluster_col]])
  cell_types  <- cluster_annotations[cluster_ids]
  
  ## Build metadata data frame
  meta_df <- data.frame(
    barcode                = stripped_matched,
    celltype               = cell_types,
    cluster                = cluster_ids,
    is_empty               = FALSE,
    percent_mt             = sample_meta_matched$percent.mt,
    nFeature_RNA           = sample_meta_matched$nFeature_RNA,
    nCount_RNA             = sample_meta_matched$nCount_RNA,
    doublet_score          = sample_meta_matched$doublet_score,
    doublet_classification = sample_meta_matched$doublet_classification,
    row.names              = stripped_matched,
    stringsAsFactors       = FALSE
  )
  
  ## Write filtered counts — Matrix Market format (cells x genes)
  counts_t <- t(counts_matched)  ## cells x genes
  Matrix::writeMM(counts_t, file.path(out_dir, "counts_filtered.mtx"))
  write.table(
    stripped_matched,
    file      = file.path(out_dir, "barcodes_filtered.tsv"),
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
  write.table(
    rownames(counts_matched),
    file      = file.path(out_dir, "genes.tsv"),
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
  write.csv(meta_df, file.path(out_dir, "metadata.csv"), row.names = FALSE)
  
  message("  Saved filtered counts: ", ncol(counts_matched), " genes x ", nrow(counts_t), " cells")
  
  ## Process unfiltered h5 if available
  if (has_unfiltered) {
    
    unfiltered_path <- file.path(
      unfiltered_h5_dir,
      paste0(sample_id, unfiltered_suffix)
    )
    
    if (file.exists(unfiltered_path)) {
      message("  Reading unfiltered h5...")
      counts_unfiltered <- Read10X_h5(unfiltered_path)
      
      ## Empty droplets = barcodes in unfiltered but NOT in filtered
      all_filtered_barcodes <- colnames(counts_filtered)
      empty_barcodes        <- setdiff(colnames(counts_unfiltered), all_filtered_barcodes)
      
      message("  Empty droplets found: ", length(empty_barcodes))
      
      counts_empty <- counts_unfiltered[, empty_barcodes]
      counts_empty_t <- t(counts_empty)  ## cells x genes
      
      ## Write empty droplet counts
      Matrix::writeMM(counts_empty_t, file.path(out_dir, "counts_empty.mtx"))
      write.table(
        empty_barcodes,
        file      = file.path(out_dir, "barcodes_empty.tsv"),
        row.names = FALSE,
        col.names = FALSE,
        quote     = FALSE
      )
      message("  Saved empty droplet counts")
      
    } else {
      message("  WARNING: unfiltered h5 not found — skipping empty droplets")
    }
  }
  
  message("  Done: ", out_dir)
}

message("\nCellSweep prep complete!")
message("Next step: run cellsweep_convert.py then cellsweep_run.py")


# Session info ------------------------------------------------------------

sessionInfo()