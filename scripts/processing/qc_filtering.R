# QC Filtering ------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Loads Cell Ranger .h5 files for all 8 13LGS ileum samples
#              (4 male, 4 female), computes per-cell percent mitochondrial
#              reads using the Ensembl 113 EnsDb (SpeTri2.0), filters
#              low-quality cells, and merges into a single Seurat object.
# Input:       data/        — Cell Ranger filtered_feature_bc_matrix.h5 files
# Output:      data/processed/seurat_merged.rds     — merged filtered object
#              data/processed/metadata_all.rds       — unfiltered metadata
#              data/processed/filter_summary.rds     — per-sample filter stats
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## AnnotationHub ID for Ensembl 113 EnsDb (Ictidomys tridecemlineatus, SpeTri2.0)
ensdb_id     <- "AH119327"

## Seurat object creation thresholds
min_cells    <- 3
min_features <- 200

## QC filtering thresholds
## mt_cutoff based on Scavuzzo et al. 2023 (doi: 10.1101/2023.06.07.544052)
mt_cutoff    <- 5     # maximum percent mitochondrial reads
nfeature_min <- 200   # minimum genes per cell (empty droplet filter)
nfeature_max <- 5000  # maximum genes per cell (doublet filter)


# Libraries ---------------------------------------------------------------

library(AnnotationHub)
library(ensembldb)
library(GenomicFeatures)
library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Connect to AnnotationHub and retrieve Ensembl 113 EnsDb for 13LGS
ah  <- AnnotationHub()
edb <- ah[[ensdb_id]]

## Confirm genome build and organism
organism(edb)
genome(edb)

## Pull MT gene names — done once, reused for all samples
## 13LGS MT genes use bare gene symbols (ND1, COX1, etc.) with no mt- prefix,
## so we cannot use pattern = "^mt-" as we would for human/mouse
mt_genes      <- genes(edb, filter = SeqNameFilter("MT"))

## Keep only named genes — drops unannotated tRNAs/rRNAs with empty gene_name
## tRNAs/rRNAs are not polyadenylated and are largely not captured by 10x
mt_gene_names <- mt_genes$gene_name[mt_genes$gene_name != ""]

## Discover all .h5 files in the data/ directory automatically
h5_files <- list.files(
  path       = here("data", "raw"),
  pattern    = "\\.h5$",
  full.names = TRUE,
)

## Parse sample IDs and sex from filenames
## e.g. "13LGS_sIBA_F1_1_filtered_feature_bc_matrix.h5" -> "13LGS_sIBA_F1_1"
sample_ids <- gsub("_filtered_feature_bc_matrix\\.h5$", "", basename(h5_files))


# Analysis ----------------------------------------------------------------

## Initialize empty lists for per-sample objects and metadata
seurat_list   <- list()
metadata_list <- list()

for (i in seq_along(h5_files)) {
  
  sample_id <- sample_ids[i]
  message("Processing: ", sample_id)
  
  ## Load Cell Ranger filtered feature-barcode matrix
  counts <- Read10X_h5(h5_files[i])
  
  ## Build Seurat object
  seurat_obj <- CreateSeuratObject(
    counts       = counts,
    project      = sample_id,
    min.cells    = min_cells,
    min.features = min_features
  )
  
  ## Add sex metadata parsed from sample ID
  seurat_obj$sex <- ifelse(grepl("_F", sample_id), "Female", "Male")
  
  ## Compute percent mitochondrial reads per cell
  ## features = used instead of pattern = because 13LGS MT genes lack mt- prefix
  seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    features = mt_gene_names[mt_gene_names %in% rownames(seurat_obj)],
    col.name = "percent.mt"
  )
  
  ## Save unfiltered metadata BEFORE filtering for QC plots
  metadata_list[[sample_id]] <- seurat_obj@meta.data |>
    transform(
      sample_id = sample_id,
      sex       = ifelse(grepl("_F", sample_id), "Female", "Male")
    )
  
  ## Filter low-quality cells
  ## - percent.mt < mt_cutoff      : removes damaged/dying cells
  ## - nFeature_RNA > nfeature_min : removes empty droplets
  ## - nFeature_RNA < nfeature_max : removes likely doublets
  seurat_obj <- subset(
    seurat_obj,
    subset = percent.mt   <  mt_cutoff    &
      nFeature_RNA >  nfeature_min &
      nFeature_RNA <  nfeature_max
  )
  
  seurat_list[[sample_id]] <- seurat_obj
  
  message(
    sample_id, ": ",
    nrow(metadata_list[[sample_id]]), " cells before filtering, ",
    ncol(seurat_obj), " after"
  )
  
}

## Merge all filtered Seurat objects into one combined object
## add.cell.ids prefixes barcodes with sample ID to avoid collisions
seurat_merged <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = sample_ids
)

## Combine unfiltered metadata for QC plotting
metadata_all <- do.call(rbind, metadata_list)

## Build per-sample filtering summary
filter_summary <- do.call(rbind, lapply(sample_ids, function(sid) {
  n_before  <- nrow(metadata_list[[sid]])
  n_after   <- ncol(seurat_list[[sid]])
  n_removed <- n_before - n_after
  data.frame(
    sample_id   = sid,
    sex         = ifelse(grepl("_F", sid), "Female", "Male"),
    Kept        = n_after,
    Removed     = n_removed,
    pct_removed = round(n_removed / n_before * 100, 1)
  )
}))


# Save outputs ------------------------------------------------------------

saveRDS(seurat_merged,  file = here("data", "processed", "seurat_merged.rds"))
saveRDS(metadata_all,   file = here("data", "processed", "metadata_all.rds"))
saveRDS(filter_summary, file = here("data", "processed", "filter_summary.rds"))


# Session info ------------------------------------------------------------

sessionInfo()