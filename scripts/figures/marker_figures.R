# Marker Gene Figures -----------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Generates marker gene visualizations including a dot plot of
#              known ileum cell type markers across clusters and a heatmap
#              of the top 10 differentially expressed genes per cluster.
# Input:       data/processed/seurat_processed.rds
#              data/processed/top10_markers.rds
# Output:      plots/dotplot_ileum_markers.pdf
#              plots/heatmap_top10_markers.pdf
# Note:        This project uses renv for reproducibility.
#              Run renv::restore() before executing this script.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Known ileum cell type marker genes
## Based on small intestinal cell type markers from mouse/human literature
## and enteric glia markers from Scavuzzo et al. 2023
ileum_markers <- c(
  ## Enterocytes
  "FABP1", "FABP2", "APOA1", "SLC5A1",
  ## Goblet cells
  "MUC2", "TFF3", "CLCA1",
  ## Paneth cells
  "DEFA5", "DEFA6", "LYZ",
  ## Enteroendocrine cells
  "CHGA", "CHGB", "SCG2",
  ## Tuft cells
  "DCLK1", "POU2F3", "TRPM5",
  ## Intestinal stem cells
  "LGR5", "OLFM4", "ASCL2",
  ## Enteric glia (Scavuzzo et al. 2023)
  "SOX10", "GFAP", "PLP1", "S100B",
  ## Enteric neurons
  "ELAVL4", "TUBB3",
  ## Smooth muscle
  "ACTA2", "MYH11",
  ## Endothelial
  "PECAM1", "CDH5",
  ## Macrophages/immune
  "CD68", "CSF1R", "PTPRC",
  ## Fibroblasts/mesenchyme
  "VIM", "COL1A1", "PDGFRA"
)


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(here)
library(Seurat)


# Load data ---------------------------------------------------------------

## Processed Seurat object from normalize_cluster.R
seurat_processed <- readRDS(here("data", "processed", "seurat_processed.rds"))

## Top 10 marker genes per cluster from marker_genes.R
top10_markers <- readRDS(here("data", "processed", "top10_markers.rds"))


# Visualization -----------------------------------------------------------

## Dot plot of known ileum cell type markers across clusters
## Size = fraction of cells expressing gene, color = average expression
p_dotplot <- DotPlot(
  seurat_processed,
  features = ileum_markers,
  assay    = "SCT"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Known ileum marker gene expression by cluster")

## Heatmap of top 10 DE genes per cluster
p_heatmap <- DoHeatmap(
  seurat_processed,
  features = top10_markers$gene,
  assay    = "SCT"
) + labs(title = "Top 10 marker genes per cluster")


# Save outputs ------------------------------------------------------------

ggsave(
  filename = here("plots", "dotplot_ileum_markers.pdf"),
  plot     = p_dotplot,
  width    = 16,
  height   = 8
)
ggsave(
  filename = here("plots", "heatmap_top10_markers.pdf"),
  plot     = p_heatmap,
  width    = 14,
  height   = 10
)


# Session info ------------------------------------------------------------

sessionInfo()