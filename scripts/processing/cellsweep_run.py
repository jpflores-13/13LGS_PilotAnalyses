#!/usr/bin/env python3
"""
CellSweep Run
-------------
Author:      JP Flores
Date:        2026-05-08
Project:     13LGS_PilotAnalyses
Description: Runs cellsweep denoising on per-sample AnnData .h5ad files
             produced by cellsweep_convert.py. Designed to be called by a
             SLURM job array, one job per sample. After denoising, subsets
             the output to non-empty barcodes only before saving, since
             empty droplets are not needed downstream.
Input:       data/processed/cellsweep/<sample_id>_input.h5ad
Output:      data/processed/cellsweep/<sample_id>_denoised.h5ad
"""

import os
import sys
import glob
import anndata as ad
import cellsweep


# Parameters ---------------------------------------------------------------

input_dir  = "data/processed/cellsweep"
output_dir = "data/processed/cellsweep"


# Identify sample from SLURM array task ID --------------------------------

input_files = sorted(glob.glob(os.path.join(input_dir, "*_input.h5ad")))

if len(input_files) == 0:
    print(f"ERROR: No _input.h5ad files found in {input_dir}")
    sys.exit(1)

task_id    = int(os.environ.get("SLURM_ARRAY_TASK_ID", 1)) - 1
input_path = input_files[task_id]
sample_id  = os.path.basename(input_path).replace("_input.h5ad", "")

print(f"Denoising sample {task_id + 1} of {len(input_files)}: {sample_id}")

n_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
print(f"Using {n_threads} thread(s)")


# Run cellsweep -----------------------------------------------------------

## Write a temporary full output (cells + empties) then subset
tmp_path = os.path.join(output_dir, f"{sample_id}_denoised_tmp.h5ad")
out_path = os.path.join(output_dir, f"{sample_id}_denoised.h5ad")
log_path = os.path.join(output_dir, f"{sample_id}_cellsweep.log")

cellsweep.denoise_count_matrix(
    adata     = input_path,
    adata_out = tmp_path,
    round_X   = True,
    threads   = n_threads,
    log_file  = log_path,
)

## Subset to real cells only — empties not needed downstream
adata_denoised = ad.read_h5ad(tmp_path)
adata_cells    = adata_denoised[~adata_denoised.obs["is_empty"]].copy()

adata_cells.write_h5ad(out_path)
os.remove(tmp_path)

print(f"\nSaved: {out_path}")
print(f"  Denoised cells: {adata_cells.shape[0]} x {adata_cells.shape[1]}")
