# QC Plots ----------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-08
# Project:     13LGS_PilotAnalyses
# Description: Generates QC plots from pre- and post-filtering metadata.
#              Produces a faceted percent.mt histogram, a stacked bar plot
#              showing cells retained vs removed, and violin plots for key
#              QC metrics (nFeature_RNA, nCount_RNA, percent.mt, doublet
#              score) per sample colored by sex.
# Input:       data/processed/metadata_all.rds
#              data/processed/filter_summary.rds
#              data/processed/metadata_processed.rds  — for post-filter violins
# Output:      plots/qc_mito_histogram.pdf
#              plots/qc_filter_summary.pdf
#              plots/qc_violin_prefilter.pdf
#              plots/qc_violin_postfilter.pdf
#              plots/qc_violin_doublets.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## MT cutoff used in qc_filtering.R
mt_cutoff    <- 5
min_features <- 200
max_features <- 5000


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(here)


# Load data ---------------------------------------------------------------

metadata_all       <- readRDS(here("data", "processed", "metadata_all.rds"))
filter_summary     <- readRDS(here("data", "processed", "filter_summary.rds"))
metadata_processed <- readRDS(here("data", "processed", "metadata_processed.rds"))


# Shared theme ------------------------------------------------------------

qc_theme <- theme_bw(base_size = 9) +
  theme(
    plot.title      = element_blank(),
    plot.caption    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    panel.grid      = element_line(color = "grey95")
  )

sex_colors <- c("Female" = "hotpink3", "Male" = "steelblue")


# Visualization — MT histogram --------------------------------------------

p_hist <- metadata_all |>
  ggplot(aes(x = percent.mt, fill = sex)) +
  geom_histogram(binwidth = 0.5, color = "grey40") +
  scale_fill_manual(values = sex_colors) +
  geom_vline(
    xintercept = mt_cutoff,
    color      = "firebrick",
    linetype   = "dashed",
    linewidth  = 0.8
  ) +
  annotate(
    geom  = "text",
    x     = mt_cutoff + 0.3,
    y     = Inf,
    label = paste0(mt_cutoff, "% cutoff"),
    color = "firebrick",
    hjust = 0,
    vjust = 2,
    size  = 2.5
  ) +
  facet_wrap(~ sample_id, ncol = 4) +
  labs(x = "% Mitochondrial Reads", y = "# of Cells", fill = "Sex") +
  qc_theme


# Visualization — filter summary bar plot ---------------------------------

filter_summary_long <- reshape(
  filter_summary,
  varying   = c("Kept", "Removed"),
  v.names   = "n_cells",
  timevar   = "status",
  times     = c("Kept", "Removed"),
  direction = "long"
)

p_filter <- filter_summary_long |>
  ggplot(aes(x = sample_id, y = n_cells, fill = interaction(status, sex))) +
  geom_col() +
  geom_text(
    data = filter_summary,
    aes(
      x     = sample_id,
      y     = Kept + Removed,
      label = paste0(pct_removed, "% removed"),
      fill  = NULL
    ),
    vjust = -0.5,
    size  = 2.5,
    color = "grey30"
  ) +
  scale_fill_manual(
    values = c(
      "Kept.Female"    = "hotpink3",
      "Removed.Female" = "pink1",
      "Kept.Male"      = "steelblue",
      "Removed.Male"   = "lightblue"
    ),
    labels = c(
      "Kept.Female"    = "Female - Kept",
      "Removed.Female" = "Female - Removed",
      "Kept.Male"      = "Male - Kept",
      "Removed.Male"   = "Male - Removed"
    )
  ) +
  labs(x = "Sample", y = "# of Cells", fill = NULL) +
  qc_theme


# Visualization — pre-filter violin plots ---------------------------------

## nFeature_RNA — genes detected per cell
p_vln_features <- metadata_all |>
  ggplot(aes(x = sample_id, y = nFeature_RNA, fill = sex)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_hline(
    yintercept = c(min_features, max_features),
    color      = "firebrick",
    linetype   = "dashed",
    linewidth  = 0.5
  ) +
  scale_fill_manual(values = sex_colors) +
  labs(x = NULL, y = "Genes per cell", fill = "Sex") +
  qc_theme

## nCount_RNA — UMIs per cell
p_vln_counts <- metadata_all |>
  ggplot(aes(x = sample_id, y = nCount_RNA, fill = sex)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  scale_fill_manual(values = sex_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "UMIs per cell (log10)", fill = "Sex") +
  qc_theme

## percent.mt — mitochondrial reads
p_vln_mt <- metadata_all |>
  ggplot(aes(x = sample_id, y = percent.mt, fill = sex)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_hline(
    yintercept = mt_cutoff,
    color      = "firebrick",
    linetype   = "dashed",
    linewidth  = 0.5
  ) +
  scale_fill_manual(values = sex_colors) +
  labs(x = NULL, y = "% Mitochondrial Reads", fill = "Sex") +
  qc_theme

## Combine pre-filter violins with patchwork
library(patchwork)
p_vln_prefilter <- (p_vln_features / p_vln_counts / p_vln_mt) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# Visualization — post-filter violin plots --------------------------------

## nFeature_RNA after filtering
p_vln_features_post <- metadata_processed |>
  ggplot(aes(x = orig.ident, y = nFeature_RNA, fill = sex)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  scale_fill_manual(values = sex_colors) +
  labs(x = NULL, y = "Genes per cell", fill = "Sex") +
  qc_theme

## nCount_RNA after filtering
p_vln_counts_post <- metadata_processed |>
  ggplot(aes(x = orig.ident, y = nCount_RNA, fill = sex)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  scale_fill_manual(values = sex_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "UMIs per cell (log10)", fill = "Sex") +
  qc_theme

## percent.mt after filtering
p_vln_mt_post <- metadata_processed |>
  ggplot(aes(x = orig.ident, y = percent.mt, fill = sex)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  scale_fill_manual(values = sex_colors) +
  labs(x = NULL, y = "% Mitochondrial Reads", fill = "Sex") +
  qc_theme

## Combine post-filter violins
p_vln_postfilter <- (p_vln_features_post / p_vln_counts_post / p_vln_mt_post) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# Visualization — doublet score violin ------------------------------------

## Only plot if doublet_score column exists in metadata
if ("doublet_score" %in% colnames(metadata_processed)) {
  
  p_vln_doublets <- metadata_processed |>
    ggplot(aes(x = orig.ident, y = doublet_score, fill = sex)) +
    geom_violin(scale = "width", linewidth = 0.3) +
    scale_fill_manual(values = sex_colors) +
    labs(x = "Sample", y = "Doublet Score (scDblFinder)", fill = "Sex") +
    qc_theme
  
  ggsave(
    filename = here("plots", "qc_violin_doublets.pdf"),
    plot     = p_vln_doublets,
    width    = 10,
    height   = 4
  )
  message("Doublet score violin saved")
  
} else {
  message("doublet_score not found in metadata — skipping doublet violin")
}


# Save outputs ------------------------------------------------------------

ggsave(
  filename = here("plots", "qc_mito_histogram.pdf"),
  plot     = p_hist,
  width    = 14,
  height   = 6
)
ggsave(
  filename = here("plots", "qc_filter_summary.pdf"),
  plot     = p_filter,
  width    = 10,
  height   = 6
)
ggsave(
  filename = here("plots", "qc_violin_prefilter.pdf"),
  plot     = p_vln_prefilter,
  width    = 12,
  height   = 10
)
ggsave(
  filename = here("plots", "qc_violin_postfilter.pdf"),
  plot     = p_vln_postfilter,
  width    = 12,
  height   = 10
)

message("QC plots saved to plots/")


# Session info ------------------------------------------------------------

sessionInfo()