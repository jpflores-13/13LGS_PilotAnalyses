#!/usr/bin/env python3
"""
CellSweep Convert
-----------------
Author:      JP Flores
Date:        2026-05-08
Project:     13LGS_PilotAnalyses
Description: Converts per-sample sparse count matrices exported by
             cellsweep_prep.R into AnnData .h5ad files ready for cellsweep
             denoising. Each .h5ad contains both filtered cells (with cell
             type labels) and empty droplets (celltype = "empty",
             is_empty = True), matching the structure expected by
             cellsweep.denoise_count_matrix().
Input:       data/processed/cellsweep/<sample_id>/
               counts.mtx
               barcodes.tsv
               genes.tsv
               metadata.csv
Output:      data/processed/cellsweep/<sample_id>_input.h5ad
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad


# Parameters ---------------------------------------------------------------

input_base_dir = "data/processed/cellsweep"
output_dir     = "data/processed/cellsweep"


# Identify sample from SLURM array task ID --------------------------------

sample_dirs = sorted([
    d for d in glob.glob(os.path.join(input_base_dir, "*"))
    if os.path.isdir(d)
])

if len(sample_dirs) == 0:
    print(f"ERROR: No sample directories found in {input_base_dir}")
    sys.exit(1)

task_id    = int(os.environ.get("SLURM_ARRAY_TASK_ID", 1)) - 1
sample_dir = sample_dirs[task_id]
sample_id  = os.path.basename(sample_dir)

print(f"Converting sample {task_id + 1} of {len(sample_dirs)}: {sample_id}")


# Load flat files ---------------------------------------------------------

print("Loading c
