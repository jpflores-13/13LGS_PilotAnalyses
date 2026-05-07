# Normalization and Clustering --------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Takes the merged filtered Seurat object, runs SCTransform
#              normalization, PCA with automatic PC selection, Harmony batch
#              correction (by sample and sex), UMAP, and clustering at
#              multiple resolutions.
# Input:       data/processed/seurat_merged.rds
# Output:      data/processed/seurat_processed.rds — normalized, clustered object
#              data/processed/elbow_df.rds          — PC selection data frame
# Note:        This project uses renv for reproducibility.
#              Run renv::restore() before executing this script.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Variables to regress out during SCTransform
vars_to_regress <- "percent.mt"

## Harmony batch correction grouping variables
harmony_vars <- c("orig.ident", "sex")

## Clustering resolutions to test
cluster_resolutions <- c(0.3, 0.5, 0.8)

## Default resolution for downstream analysis (middle of the three)
default_resolution <- cluster_resolutions[2]


# Libraries ---------------------------------------------------------------

library(harmony)
library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Merged filtered Seurat object from qc_filtering.R
seurat_merged <- readRDS(here("data", "processed", "seurat_merged.rds"))


# Analysis ----------------------------------------------------------------

## Join layers before SCTransform (required for Seurat v5 merged objects)
seurat_merged <- JoinLayers(seurat_merged)

## SCTransform normalization — regresses out percent.mt
## return.only.var.genes = FALSE retains all genes for downstream flexibility
## Increase future globals size limit for SCTransform on large objects
## Default 500 MiB is too small for merged multi-sample Seurat objects
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GiB
seurat_merged <- SCTransform(
  seurat_merged,
  vars.to.regress       = vars_to_regress,
  return.only.var.genes = FALSE,
  verbose               = FALSE
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
  pct  = pct_stdv,
  cumu = cumulative,
  rank = seq_along(pct_stdv),
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