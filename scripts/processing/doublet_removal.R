# Doublet Removal ---------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-08
# Project:     13LGS_PilotAnalyses
# Description: Runs scDblFinder on a single sample to identify and flag
#              putative doublets. Designed to be called by a SLURM job array,
#              one job per sample. Saves two outputs per sample:
#                1. Flagged object — all cells with doublet scores in metadata
#                2. Cleaned object — singlets only, ready for SCTransform
#              scDblFinder is used instead of DoubletFinder for better HPC
#              compatibility and simpler API.
# Input:       data/processed/seurat_merged.rds
# Output:      data/processed/doublets/<sample_id>_flagged.rds
#              data/processed/doublets/<sample_id>_clean.rds
# Note:        Install scDblFinder with: BiocManager::install("scDblFinder")
#              Submit as a SLURM job array, not interactively.This project uses renv for reproducibility.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Random seed for reproducibility
set_seed <- 42 #Jackie Robinson!!!


# Libraries ---------------------------------------------------------------

library(here)
library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)


# Load data ---------------------------------------------------------------

## Get SLURM array task ID to determine which sample to process
task_id   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))

## Load merged object to extract sample IDs
seurat_merged <- readRDS(here("data", "processed", "seurat_merged.rds"))
sample_ids    <- unique(seurat_merged$orig.ident)
sample_id     <- sample_ids[task_id]

message("Processing doublets for sample ", task_id, " of ", length(sample_ids), ": ", sample_id)

## Subset to this sample
seurat_obj <- subset(seurat_merged, subset = orig.ident == sample_id)
n_cells    <- ncol(seurat_obj)

message("  Cells before doublet removal: ", n_cells)


# Analysis ----------------------------------------------------------------

## Convert Seurat object to SingleCellExperiment for scDblFinder
## scDblFinder works natively with SCE objects
sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")

## Run scDblFinder
## scDblFinder automatically estimates doublet rate from cell count
set.seed(set_seed)
sce <- scDblFinder(sce)

## Extract doublet calls and scores back into Seurat metadata
seurat_obj$doublet_classification <- ifelse(
  sce$scDblFinder.class == "doublet", "Doublet", "Singlet"
)
seurat_obj$doublet_score <- sce$scDblFinder.score

## Summarize results
n_doublets  <- sum(seurat_obj$doublet_classification == "Doublet")
n_singlets  <- sum(seurat_obj$doublet_classification == "Singlet")
pct_doublet <- round(n_doublets / n_cells * 100, 1)

message("  Doublets found:    ", n_doublets, " (", pct_doublet, "%)")
message("  Singlets retained: ", n_singlets)

## Flagged object — all cells with doublet classification in metadata
seurat_flagged <- seurat_obj

## Cleaned object — singlets only, ready for SCTransform
seurat_clean <- subset(seurat_obj, subset = doublet_classification == "Singlet")


# Save outputs ------------------------------------------------------------

## Create output directory if it doesn't exist
dir.create(here("data", "processed", "doublets"), showWarnings = FALSE)

saveRDS(
  seurat_flagged,
  file = here("data", "processed", "doublets", paste0(sample_id, "_flagged.rds"))
)
saveRDS(
  seurat_clean,
  file = here("data", "processed", "doublets", paste0(sample_id, "_clean.rds"))
)

message("Saved:")
message("  ", sample_id, "_flagged.rds — all cells with doublet scores")
message("  ", sample_id, "_clean.rds   — singlets only (", n_singlets, " cells)")


# Session info ------------------------------------------------------------

sessionInfo()