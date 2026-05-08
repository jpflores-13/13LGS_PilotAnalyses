# Marker Plots ------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Generates marker gene plots inspired by Scavuzzo et al. 2026
#              style. Dot plot has genes on y-axis and clusters on x-axis
#              with grey shading boxes grouping cell types. Heatmap shows
#              top 10 marker genes per cluster with colored annotation bar
#              on top. Volcano plot annotates top genes with text labels.
#              Set input_suffix to "_denoised" for cellsweep pass plots.
# Input:       data/processed/marker_expr<input_suffix>.rds
#              data/processed/metadata_processed<input_suffix>.rds
#              data/processed/seurat_markers<input_suffix>.rds
#              data/processed/top10_markers<input_suffix>.rds
# Output:      plots/dotplot_ileum_markers<input_suffix>.pdf
#              plots/heatmap_top10_markers<input_suffix>.pdf
#              plots/volcano_markers<input_suffix>.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## Input/output suffix — set to "_denoised" for cellsweep pass
input_suffix <- "_denoised"

## Default clustering resolution
default_resolution <- 0.5
cluster_col        <- paste0("SCT_snn_res.", default_resolution)

## Number of top genes to label on volcano
n_label <- 5

## First-pass cluster annotations — must match umap_plots.R
cluster_annotations <- c(
  "0"  = "Enterocyte",
  "1"  = "Enteric Neuron",
  "2"  = "Unknown",
  "3"  = "T cell",
  "4"  = "Enterocyte",
  "5"  = "Fibroblast",
  "6"  = "Enteric Neuron",
  "7"  = "Immune",
  "8"  = "Enterocyte",
  "9"  = "Goblet",
  "10" = "T cell",
  "11" = "Smooth Muscle",
  "12" = "Endothelial",
  "13" = "Lymphatic Endothelial",
  "14" = "Enteric Neuron",
  "15" = "Fibroblast",
  "16" = "Enteric Neuron",
  "17" = "Enteric Neuron",
  "18" = "Enteric Neuron",
  "19" = "Enteric Glia",
  "20" = "Enteric Neuron",
  "21" = "Enteric Neuron",
  "22" = "Enteric Neuron",
  "23" = "Enteric Glia",
  "24" = "Smooth Muscle",
  "25" = "Smooth Muscle",
  "26" = "Smooth Muscle",
  "27" = "Unknown",
  "28" = "Smooth Muscle"
)

## Cell type color palette
celltype_colors <- c(
  "Enterocyte"            = "#4E9AF1",
  "Goblet"                = "#6DC8A0",
  "T cell"                = "#06D6A0",
  "Immune"                = "#2EC4B6",
  "Enteric Neuron"        = "#C77DFF",
  "Enteric Glia"          = "#F15BB5",
  "Smooth Muscle"         = "#FF6B6B",
  "Endothelial"           = "#FFD166",
  "Lymphatic Endothelial" = "#F4A261",
  "Fibroblast"            = "#118AB2",
  "Unknown"               = "grey70"
)

## Dot plot gene groups and cell type assignments
gene_groups <- data.frame(
  gene  = c(
    "FABP1", "FABP2", "APOA1", "SLC5A1",
    "MUC2", "TFF3", "CLCA1",
    "DEFA5", "DEFA6", "LYZ",
    "CHGA", "CHGB", "SCG2",
    "DCLK1", "POU2F3", "TRPM5",
    "LGR5", "OLFM4", "ASCL2",
    "SOX10", "GFAP", "PLP1", "S100B",
    "ELAVL4", "TUBB3",
    "ACTA2", "MYH11",
    "PECAM1", "CDH5",
    "CD68", "CSF1R", "PTPRC",
    "VIM", "COL1A1", "PDGFRA"
  ),
  group = c(
    rep("Enterocyte",      4),
    rep("Goblet",          3),
    rep("Paneth",          3),
    rep("Enteroendocrine", 3),
    rep("Tuft",            3),
    rep("ISC",             3),
    rep("Enteric Glia",    4),
    rep("Enteric Neuron",  2),
    rep("Smooth Muscle",   2),
    rep("Endothelial",     2),
    rep("Immune",          3),
    rep("Fibroblast",      3)
  ),
  stringsAsFactors = FALSE
)

## Group display order and colors for shading boxes
group_levels <- c(
  "Enterocyte", "Goblet", "Paneth", "Enteroendocrine", "Tuft", "ISC",
  "Enteric Glia", "Enteric Neuron", "Smooth Muscle", "Endothelial",
  "Immune", "Fibroblast"
)

group_colors <- c(
  "Enterocyte"      = "#4E9AF1",
  "Goblet"          = "#6DC8A0",
  "Paneth"          = "#A8D8A8",
  "Enteroendocrine" = "#F4A261",
  "Tuft"            = "#E76F51",
  "ISC"             = "#9B5DE5",
  "Enteric Glia"    = "#F15BB5",
  "Enteric Neuron"  = "#C77DFF",
  "Smooth Muscle"   = "#FF6B6B",
  "Endothelial"     = "#FFD166",
  "Immune"          = "#06D6A0",
  "Fibroblast"      = "#118AB2"
)


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(here)


# Load data ---------------------------------------------------------------

marker_expr        <- readRDS(here("data", "processed", paste0("marker_expr",        input_suffix, ".rds")))
metadata_processed <- readRDS(here("data", "processed", paste0("metadata_processed", input_suffix, ".rds")))
seurat_markers     <- readRDS(here("data", "processed", paste0("seurat_markers",     input_suffix, ".rds")))
top10_markers      <- readRDS(here("data", "processed", paste0("top10_markers",      input_suffix, ".rds")))

cluster_labels <- metadata_processed[[cluster_col]]
cluster_order  <- as.character(sort(as.numeric(unique(cluster_labels))))


# Wrangle data — dot plot -------------------------------------------------

markers_present <- gene_groups$gene[gene_groups$gene %in% colnames(marker_expr)]

if (length(markers_present) < nrow(gene_groups)) {
  message(
    nrow(gene_groups) - length(markers_present),
    " marker genes not detected: ",
    paste(setdiff(gene_groups$gene, markers_present), collapse = ", ")
  )
}

## Build dot plot summary per cluster
dot_df <- do.call(rbind, lapply(cluster_order, function(cl) {
  cells <- as.character(cluster_labels) == cl
  do.call(rbind, lapply(markers_present, function(gene) {
    expr <- marker_expr[cells, gene]
    data.frame(
      cluster     = cl,
      gene        = gene,
      avg_expr    = mean(expr),
      pct_express = mean(expr > 0) * 100,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df         <- merge(dot_df, gene_groups, by = "gene")
gene_order     <- gene_groups$gene[gene_groups$gene %in% markers_present]
dot_df$gene    <- factor(dot_df$gene,    levels = rev(gene_order))
dot_df$cluster <- factor(dot_df$cluster, levels = cluster_order)
dot_df$group   <- factor(dot_df$group,   levels = group_levels)

## Build grey shading rectangles for gene groups on dot plot
group_positions <- do.call(rbind, lapply(group_levels, function(g) {
  genes_in_group <- rev(gene_order)[rev(gene_order) %in%
                                      gene_groups$gene[gene_groups$group == g]]
  if (length(genes_in_group) == 0) return(NULL)
  positions <- which(levels(dot_df$gene) %in% genes_in_group)
  data.frame(
    group  = g,
    ymin   = min(positions) - 0.5,
    ymax   = max(positions) + 0.5,
    stringsAsFactors = FALSE
  )
}))
group_positions$fill <- ifelse(
  seq_len(nrow(group_positions)) %% 2 == 0, "grey92", "white"
)

## Gene label colors by group
gene_color_map <- setNames(
  group_colors[gene_groups$group[match(gene_order, gene_groups$gene)]],
  gene_order
)


# Wrangle data — heatmap --------------------------------------------------

heatmap_genes <- seurat_markers$gene[seurat_markers$gene %in% colnames(marker_expr)]
heatmap_genes <- unique(heatmap_genes)

message("Building heatmap with ", length(heatmap_genes), " genes (all significant markers)")

## Average expression per cluster per gene
heatmap_df <- do.call(rbind, lapply(cluster_order, function(cl) {
  cells <- as.character(cluster_labels) == cl
  do.call(rbind, lapply(heatmap_genes, function(gene) {
    data.frame(
      cluster  = cl,
      gene     = gene,
      avg_expr = mean(marker_expr[cells, gene]),
      stringsAsFactors = FALSE
    )
  }))
}))

## Reshape to wide matrix for z-scoring
heatmap_wide <- reshape(
  heatmap_df,
  idvar     = "gene",
  timevar   = "cluster",
  direction = "wide"
)
rownames(heatmap_wide) <- heatmap_wide$gene
heatmap_wide$gene      <- NULL
colnames(heatmap_wide) <- gsub("avg_expr\\.", "", colnames(heatmap_wide))
heatmap_wide           <- heatmap_wide[, cluster_order]

## Z-score each gene across clusters
heatmap_mat <- as.matrix(heatmap_wide)
heatmap_z   <- t(scale(t(heatmap_mat)))

## Order genes by cluster of peak expression for diagonal pattern
peak_cluster    <- apply(heatmap_z, 1, which.max)
gene_order_heat <- rownames(heatmap_z)[order(peak_cluster)]

## Reshape back to long for ggplot
heatmap_z_long <- data.frame(
  gene    = rep(rownames(heatmap_z), times = ncol(heatmap_z)),
  cluster = rep(cluster_order, each = nrow(heatmap_z)),
  z_score = as.vector(heatmap_z[, cluster_order])
)

heatmap_z_long$gene    <- factor(heatmap_z_long$gene,    levels = gene_order_heat)
heatmap_z_long$cluster <- factor(heatmap_z_long$cluster, levels = cluster_order)

## Build colored annotation bar on top of heatmap
anno_bar <- data.frame(
  cluster   = cluster_order,
  cell_type = cluster_annotations[cluster_order],
  stringsAsFactors = FALSE
)
anno_bar$cluster   <- factor(anno_bar$cluster,   levels = cluster_order)
anno_bar$cell_type <- factor(anno_bar$cell_type, levels = names(celltype_colors))
anno_bar$y         <- 1


# Visualization — dot plot ------------------------------------------------

p_dotplot <- ggplot() +
  geom_rect(
    data = group_positions,
    aes(ymin = ymin, ymax = ymax, fill = fill),
    xmin = -Inf, xmax = Inf,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  geom_point(
    data = dot_df,
    aes(x = cluster, y = gene, size = pct_express, color = avg_expr)
  ) +
  scale_color_viridis_c(option = "viridis", name = "avg exp") +
  scale_size_continuous(
    range  = c(0.2, 5),
    breaks = c(0, 25, 50, 75),
    name   = "% exp"
  ) +
  geom_text(
    data = group_positions |>
      transform(
        label = group,
        x     = as.numeric(factor(cluster_order[length(cluster_order)],
                                  levels = cluster_order)) + 0.8,
        y     = (ymin + ymax) / 2
      ),
    aes(x = length(cluster_order) + 0.8, y = (ymin + ymax) / 2, label = group),
    size        = 2.5,
    hjust       = 0,
    color       = group_colors[group_positions$group],
    fontface    = "bold",
    inherit.aes = FALSE
  ) +
  labs(x = "Cluster", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 9) +
  theme(
    plot.title      = element_blank(),
    plot.caption    = element_blank(),
    axis.text.y     = element_text(
      size  = 7,
      color = gene_color_map[rev(gene_order)]
    ),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid      = element_line(color = "grey97"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(5, 80, 5, 5)
  )


# Visualization — heatmap -------------------------------------------------

p_anno <- anno_bar |>
  ggplot(aes(x = cluster, y = y, fill = cell_type)) +
  geom_tile() +
  scale_fill_manual(values = celltype_colors, name = "Cell Type") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text     = element_text(size = 7),
    legend.title    = element_text(size = 8)
  )

p_heatmap_main <- heatmap_z_long |>
  ggplot(aes(x = cluster, y = gene, fill = z_score)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(
    colors = c("#440154", "#31688e", "#35b779", "#fde725"),
    name   = "Z-score",
    limits = c(-3, 3),
    oob    = scales::squish
  ) +
  labs(x = "Cluster", y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    plot.title   = element_blank(),
    plot.caption = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    legend.key.height = unit(0.8, "cm")
  )

library(patchwork)
p_heatmap <- p_anno / p_heatmap_main +
  plot_layout(heights = c(0.03, 1), guides = "collect")


# Visualization — volcano -------------------------------------------------

if (nrow(seurat_markers) > 0) {
  
  top_labels <- seurat_markers |>
    (\(df) split(df, df$cluster))() |>
    lapply(\(x) x[order(-x$avg_log2FC), ][1:min(n_label, nrow(x)), ]) |>
    do.call(what = rbind)
  
  p_volcano <- seurat_markers |>
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300))) +
    geom_point(size = 0.3, alpha = 0.4, color = "grey70") +
    geom_point(
      data  = top_labels,
      aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300)),
      color = "steelblue",
      size  = 0.8
    ) +
    ggrepel::geom_text_repel(
      data         = top_labels,
      aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300), label = gene),
      size         = 1.8,
      color        = "steelblue",
      fontface     = "italic",
      box.padding  = 0.2,
      max.overlaps = 10,
      segment.size = 0.2,
      segment.color = "grey60"
    ) +
    facet_wrap(~ cluster, ncol = 5, scales = "free") +
    labs(
      x = "Average log2 Fold Change",
      y = "-log10 Adjusted P-value"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title       = element_blank(),
      plot.caption     = element_blank(),
      strip.text       = element_text(size = 7, face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      axis.text        = element_text(size = 6),
      panel.spacing    = unit(0.3, "lines")
    )
  
  ggsave(
    filename = here("plots", paste0("volcano_markers", input_suffix, ".pdf")),
    plot     = p_volcano,
    width    = 16,
    height   = 14
  )
  message("Volcano saved")
  
} else {
  message("seurat_markers empty — skipping volcano")
}


# Save outputs ------------------------------------------------------------

ggsave(
  filename = here("plots", paste0("dotplot_ileum_markers",  input_suffix, ".pdf")),
  plot     = p_dotplot,
  width    = 14,
  height   = 12
)
ggsave(
  filename = here("plots", paste0("heatmap_top10_markers", input_suffix, ".pdf")),
  plot     = p_heatmap,
  width    = 12,
  height   = 14
)

message("Marker plots saved to plots/")


# Session info ------------------------------------------------------------

sessionInfo()