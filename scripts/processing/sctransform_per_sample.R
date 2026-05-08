# SCTransform Per Sample --------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Runs SCTransform on a single sample. Designed to be called
#              by a SLURM job array, one job per sample. The SLURM_ARRAY_TASK_ID
#              environment variable is used to select which sample to process.
#              Much more memory efficient than running SCTransform on the full
#              merged object.
# Input:       data/processed/seurat_merged.rds
# Output:      data/processed/sct/<sample_id>_sct.rds — one per sample
# Note:        This project uses renv for reproducibility.
#              Run renv::restore() before executing this script.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Variables to regress out during SCTransform
vars_to_regress <- "percent.mt"

## Increase future globals size limit for SCTransform
## Default 500 MiB is too small even for single samples on some datasets
options(future.globals.maxSize = 4000 * 1024^2)  # 4 GiB


# Libraries ---------------------------------------------------------------

library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Get SLURM array task ID to determine which sample to process
## When run locally (not via SLURM), defaults to 1 for testing
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))

## Load merged object to extract sample list and IDs
seurat_merged <- readRDS(here("data", "processed", "seurat_merged.rds"))

## Get ordered sample IDs from the merged object
sample_ids <- unique(seurat_merged$orig.ident)

## Select the sample for this job
sample_id <- sample_ids[task_id]
message("Processing sample ", task_id, " of ", length(sample_ids), ": ", sample_id)

## Subset to just this sample
seurat_obj <- subset(seurat_merged, subset = orig.ident == sample_id)


# Analysis ----------------------------------------------------------------

## Run SCTransform on this single sample
## Much more memory efficient than running on the full merged object
seurat_obj <- SCTransform(
  seurat_obj,
  vars.to.regress       = vars_to_regress,
  return.only.var.genes = FALSE,
  verbose               = FALSE
)


# Save outputs ------------------------------------------------------------

## Create output directory if it doesn't exist
dir.create(here("data", "processed", "sct"), showWarnings = FALSE)

saveRDS(
  seurat_obj,
  file = here("data", "processed", "sct", paste0(sample_id, "_sct.rds"))
)

message("Done: ", sample_id)


# Session info ------------------------------------------------------------

sessionInfo()