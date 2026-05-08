# Normalization and Clustering --------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Merges per-sample SCTransform outputs from the SLURM job array,
#              runs PCA with automatic PC selection, Harmony batch correction
#              (by sample and sex), UMAP, and clustering at multiple resolutions.
#              SCTransform is run per-sample (via sctransform_per_sample.R)
#              rather than on the merged object for memory efficiency.
#              When run after cellsweep, set input_dir to sct_denoised and
#              output_suffix to "_denoised".
# Input:       data/processed/<input_dir>/<sample_id>_sct.rds — one per sample
# Output:      data/processed/seurat_processed<output_suffix>.rds
#              data/processed/elbow_df<output_suffix>.rds
# Note:        Run after all sctransform_per_sample.R jobs have completed.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Input directory — switch to "sct_denoised" for denoised counts
input_dir <- "sct_denoised"

## Output suffix — appended to output filenames to avoid overwriting first pass
output_suffix <- "_denoised"

## Harmony batch correction grouping variables
harmony_vars <- c("orig.ident", "sex")

## Clustering resolutions to test
cluster_resolutions <- c(0.3, 0.5, 0.8)

## Default resolution for downstream analysis
default_resolution <- cluster_resolutions[2]

## Increase future globals size limit for large merged objects
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GiB


# Libraries ---------------------------------------------------------------

library(harmony)
library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Discover all per-sample SCTransform outputs from the job array
sct_files <- list.files(
  path       = here("data", "processed", input_dir),
  pattern    = "_sct\\.rds$",
  full.names = TRUE
)

message("Found ", length(sct_files), " per-sample SCTransform objects")

if (length(sct_files) == 0) {
  stop(
    "No SCTransform .rds files found in data/processed/", input_dir, "/\n",
    "here() resolves to: ", here(), "\n",
    "Full expected path: ", here("data", "processed", input_dir)
  )
}

sct_list       <- lapply(sct_files, readRDS)
sct_sample_ids <- gsub("_sct\\.rds$", "", basename(sct_files))

message("Loaded objects: ", paste(sct_sample_ids, collapse = ", "))


# Analysis ----------------------------------------------------------------

## Merge all SCTransform-normalized objects into one
seurat_merged <- merge(
  x            = sct_list[[1]],
  y            = sct_list[-1],
  add.cell.ids = sct_sample_ids
)

message("Merged object: ", ncol(seurat_merged), " cells, ", nrow(seurat_merged), " features")

## Join layers after merge on RNA assay only
## SCT assay does not support JoinLayers in this per-sample SCTransform workflow
seurat_merged <- JoinLayers(seurat_merged, assay = "RNA")

## Select variable features across samples for PCA
## Uses consensus variable features across all per-sample SCTransform models
VariableFeatures(seurat_merged) <- SelectIntegrationFeatures(
  object.list = sct_list,
  nfeatures   = 3000
)

message("Selected ", length(VariableFeatures(seurat_merged)), " variable features")

## PCA on SCTransform variable features
seurat_merged <- RunPCA(
  seurat_merged,
  features = VariableFeatures(seurat_merged),
  verbose  = FALSE
)

message("PCA complete")

## Automatic PC selection (Scavuzzo et al. 2023 approach)
## co1: first PC where cumulative variance > 90% AND per-PC variance < 5%
## co2: last PC where drop in variance between adjacent PCs > 0.1%
## min_pc: minimum of co1 and co2
stdv       <- seurat_merged[["pca"]]@stdev
sum_stdv   <- sum(stdv)
pct_stdv   <- (stdv / sum_stdv) * 100
cumulative <- cumsum(pct_stdv)

co1 <- which(cumulative > 90 & pct_stdv < 5)[1]
co2 <- sort(
  which((pct_stdv[1:(length(pct_stdv) - 1)] - pct_stdv[2:length(pct_stdv)]) > 0.1),
  decreasing = TRUE
)[1] + 1

min_pc <- min(co1, co2)
message("Optimal number of PCs: ", min_pc)

## Save elbow plot data frame for figures script
elbow_df <- data.frame(
  pct    = pct_stdv,
  cumu   = cumulative,
  rank   = seq_along(pct_stdv),
  min_pc = min_pc
)

## Harmony batch correction — corrects for sample and sex effects
## Updated for Harmony2 API — assay.use replaced with reduction/reduction.save
seurat_merged <- RunHarmony(
  seurat_merged,
  group.by.vars  = harmony_vars,
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)

message("Harmony complete")

## UMAP and neighbor graph using Harmony-corrected embeddings
seurat_merged <- RunUMAP(
  seurat_merged,
  reduction = "harmony",
  dims      = 1:min_pc,
  verbose   = FALSE
)

message("UMAP complete")

seurat_merged <- FindNeighbors(
  seurat_merged,
  reduction = "harmony",
  dims      = 1:min_pc,
  verbose   = FALSE
)

## Cluster at all specified resolutions
for (res in cluster_resolutions) {
  seurat_merged <- FindClusters(
    seurat_merged,
    resolution = res,
    verbose    = FALSE
  )
  message("Clustering at resolution ", res, " complete")
}

Idents(seurat_merged) <- paste0("SCT_snn_res.", default_resolution)

message("Clustering complete — default resolution: ", default_resolution)


# Save outputs ------------------------------------------------------------

saveRDS(
  seurat_merged,
  file = here("data", "processed", paste0("seurat_processed", output_suffix, ".rds"))
)
saveRDS(
  elbow_df,
  file = here("data", "processed", paste0("elbow_df", output_suffix, ".rds"))
)

message("Saved seurat_processed", output_suffix, ".rds and elbow_df", output_suffix, ".rds")


# Session info ------------------------------------------------------------

sessionInfo()