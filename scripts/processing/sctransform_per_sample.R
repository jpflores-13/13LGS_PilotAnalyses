# SCTransform Per Sample --------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-08
# Project:     13LGS_PilotAnalyses
# Description: Runs SCTransform on a single doublet-cleaned sample. Designed
#              to be called by a SLURM job array, one job per sample. Reads
#              cleaned singlet objects from doublet_removal.R rather than the
#              raw merged object. Run after all doublet_removal.R jobs finish.
#              When run after cellsweep, set input_dir to
#              data/processed/cellsweep_clean and output_dir to
#              data/processed/sct_denoised.
# Input:       data/processed/<input_dir>/<sample_id>_clean.rds
# Output:      data/processed/<output_dir>/<sample_id>_sct.rds
# Note:        This project uses renv for reproducibility.
#              Submit as a SLURM job array after doublet_removal.R completes.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Variables to regress out during SCTransform
vars_to_regress <- "percent.mt"

## Input directory — switch to "cellsweep_clean" for denoised counts
input_dir  <- "cellsweep_clean"

## Output directory — switch to "sct_denoised" for denoised counts
output_dir <- "sct_denoised"

## Increase future globals size limit for SCTransform
## Default 500 MiB is too small even for single samples on some datasets
options(future.globals.maxSize = 4000 * 1024^2)  # 4 GiB


# Libraries ---------------------------------------------------------------

library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Get SLURM array task ID to determine which sample to process
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))

## Get ordered sample IDs from the input directory
clean_files <- list.files(
  path       = here("data", "processed", input_dir),
  pattern    = "_clean\\.rds$",
  full.names = TRUE
)
sample_ids <- gsub("_clean\\.rds$", "", basename(clean_files))
sample_id  <- sample_ids[task_id]

message("Processing SCTransform for sample ", task_id, " of ", length(sample_ids), ": ", sample_id)

## Load cleaned object
seurat_obj <- readRDS(clean_files[task_id])

message("  Cells: ", ncol(seurat_obj))


# Analysis ----------------------------------------------------------------

## Run SCTransform on this single cleaned sample
## Much more memory efficient than running on the full merged object
seurat_obj <- SCTransform(
  seurat_obj,
  vars.to.regress       = vars_to_regress,
  return.only.var.genes = FALSE,
  verbose               = FALSE
)

message("  SCTransform complete")


# Save outputs ------------------------------------------------------------

dir.create(here("data", "processed", output_dir), showWarnings = FALSE)

saveRDS(
  seurat_obj,
  file = here("data", "processed", output_dir, paste0(sample_id, "_sct.rds"))
)

message("Saved: ", sample_id, "_sct.rds to ", output_dir)


# Session info ------------------------------------------------------------

sessionInfo()