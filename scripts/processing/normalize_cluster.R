# Normalization and Clustering --------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Merges per-sample SCTransform outputs from the SLURM job array,
#              runs PCA with automatic PC selection, Harmony batch correction
#              (by sample and sex), UMAP, and clustering at multiple resolutions.
#              SCTransform is run per-sample (via sctransform_per_sample.R)
#              rather than on the merged object for memory efficiency.
# Input:       data/processed/sct/<sample_id>_sct.rds — one per sample
# Output:      data/processed/seurat_processed.rds
#              data/processed/elbow_df.rds
# Note:       
#              Run after all sctransform_per_sample.R jobs have completed.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

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
  path       = here("data", "processed", "sct"),
  pattern    = "_sct\\.rds$",
  full.names = TRUE
)

message("Found ", length(sct_files), " per-sample SCTransform objects")

## Load all per-sample objects into a list
sct_list <- lapply(sct_files, readRDS)

## Parse sample IDs from filenames for add.cell.ids
sct_sample_ids <- gsub("_sct\\.rds$", "", basename(sct_files))


# Analysis ----------------------------------------------------------------

## Merge all SCTransform-normalized objects into one
## PrepSCTFindMarkers() will handle integration of SCT models downstream
seurat_merged <- merge(
  x            = sct_list[[1]],
  y            = sct_list[-1],
  add.cell.ids = sct_sample_ids
)

## Join layers after merge (required for Seurat v5)
seurat_merged <- JoinLayers(seurat_merged)

## Select variable features across samples for PCA
## SelectIntegrationFeatures5() is the Seurat v5 way to do this post-SCT merge
VariableFeatures(seurat_merged) <- SelectIntegrationFeatures(
  object.list = sct_list,
  nfeatures   = 3000
)

## PCA on SCTransform variable features
seurat_merged <- RunPCA(
  seurat_merged,
  features = VariableFeatures(seurat_merged),
  verbose  = FALSE
)

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
seurat_merged <- RunHarmony(
  seurat_merged,
  group.by.vars = harmony_vars,
  assay.use     = "SCT",
  verbose       = FALSE
)

## UMAP and neighbor graph using Harmony-corrected embeddings
seurat_merged <- RunUMAP(
  seurat_merged,
  reduction = "harmony",
  dims      = 1:min_pc,
  verbose   = FALSE
)

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

## Set default identity to the specified default resolution
Idents(seurat_merged) <- paste0("SCT_snn_res.", default_resolution)


# Save outputs ------------------------------------------------------------

saveRDS(seurat_merged, file = here("data", "processed", "seurat_processed.rds"))
saveRDS(elbow_df,      file = here("data", "processed", "elbow_df.rds"))


# Session info ------------------------------------------------------------

sessionInfo()