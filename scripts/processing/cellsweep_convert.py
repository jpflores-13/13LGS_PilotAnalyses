#!/usr/bin/env python3
"""
CellSweep Convert
-----------------
Author:      JP Flores
Date:        2026-05-08
Project:     13LGS_PilotAnalyses
Description: Converts per-sample sparse count matrices exported by
             cellsweep_prep.R into AnnData .h5ad files for cellsweep
             denoising. Handles both filtered cells and empty droplets
             when unfiltered .h5 files were available.
Input:       data/processed/cellsweep/<sample_id>/
               counts_filtered.mtx
               barcodes_filtered.tsv
               genes.tsv
               metadata.csv
               counts_empty.mtx       (optional)
               barcodes_empty.tsv     (optional)
Output:      data/processed/cellsweep/<sample_id>_input.h5ad
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


# Parameters ---------------------------------------------------------------

input_base_dir = "data/processed/cellsweep"
output_dir     = "data/processed/cellsweep"


# Get sample from SLURM array task ID -------------------------------------

## Get all per-sample directories
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


# Load filtered cells -----------------------------------------------------

print("Loading filtered counts...")
counts_filtered = scipy.io.mmread(
  os.path.join(sample_dir, "counts_filtered.mtx")
).tocsr().astype("float32")

barcodes_filtered = pd.read_csv(
  os.path.join(sample_dir, "barcodes_filtered.tsv"),
  header = None
)[0].tolist()

genes = pd.read_csv(
  os.path.join(sample_dir, "genes.tsv"),
  header = None
)[0].tolist()

metadata = pd.read_csv(
  os.path.join(sample_dir, "metadata.csv"),
  index_col = "barcode"
)

print(f"  Filtered cells: {counts_filtered.shape[0]} x {counts_filtered.shape[1]}")

## Build filtered AnnData
adata_filtered = ad.AnnData(
  X   = counts_filtered,
  obs = metadata.loc[barcodes_filtered],
  var = pd.DataFrame(index = genes)
)
adata_filtered.obs_names = barcodes_filtered
adata_filtered.var_names = genes


# Load empty droplets if available ----------------------------------------

empty_mtx_path = os.path.join(sample_dir, "counts_empty.mtx")
empty_bc_path  = os.path.join(sample_dir, "barcodes_empty.tsv")
has_empty      = os.path.exists(empty_mtx_path) and os.path.exists(empty_bc_path)

if has_empty:
  print("Loading empty droplet counts...")
counts_empty = scipy.io.mmread(empty_mtx_path).tocsr().astype("float32")
barcodes_empty = pd.read_csv(empty_bc_path, header=None)[0].tolist()

print(f"  Empty droplets: {counts_empty.shape[0]} x {counts_empty.shape[1]}")

## Build empty droplet obs metadata
obs_empty = pd.DataFrame({
  "celltype"               : "empty",
  "cluster"                : np.nan,
  "is_empty"               : True,
  "percent_mt"             : np.nan,
  "nFeature_RNA"           : np.nan,
  "nCount_RNA"             : np.nan,
  "doublet_score"          : np.nan,
  "doublet_classification" : np.nan
}, index = barcodes_empty)

adata_empty = ad.AnnData(
  X   = counts_empty,
  obs = obs_empty,
  var = pd.DataFrame(index = genes)
)
adata_empty.obs_names = barcodes_empty
adata_empty.var_names = genes

## Concatenate filtered cells + empty droplets
print("Concatenating filtered cells + empty droplets...")
adata = ad.concat(
  [adata_filtered, adata_empty],
  axis    = 0,
  join    = "outer",
  fill_value = 0
)
adata.var_names = genes

else:
  print("No empty droplets found — using filtered cells only")
adata = adata_filtered
adata.obs["is_empty"] = False


# Write h5ad --------------------------------------------------------------

out_path = os.path.join(output_dir, f"{sample_id}_input.h5ad")
print(f"\nWriting AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  Cell types: {adata.obs['celltype'].dropna().unique().tolist()}")
print(f"  Empty droplets: {adata.obs['is_empty'].sum()}")
adata.write_h5ad(out_path)
print(f"  Saved: {out_path}")
