# QC Plots ----------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-05-07
# Project:     13LGS_PilotAnalyses
# Description: Generates QC plots from the pre-filtering metadata and
#              per-sample filtering summary. Produces a faceted percent.mt
#              histogram and a stacked bar plot showing cells retained vs
#              removed after QC filtering. Loads only lightweight .rds objects
#              so this script can be run interactively without memory issues.
# Input:       data/processed/metadata_all.rds
#              data/processed/filter_summary.rds
# Output:      plots/qc_mito_histogram.pdf
#              plots/qc_filter_summary.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

## MT cutoff used in qc_filtering.R — must match for accurate labeling
mt_cutoff <- 5


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(here)


# Load data ---------------------------------------------------------------

## Lightweight metadata objects — no Seurat object needed for these plots
metadata_all   <- readRDS(here("data", "processed", "metadata_all.rds"))
filter_summary <- readRDS(here("data", "processed", "filter_summary.rds"))


# Visualization -----------------------------------------------------------

## Faceted histogram of percent.mt per sample colored by sex
p_hist <- metadata_all |>
  ggplot(aes(x = percent.mt, fill = sex)) +
  geom_histogram(binwidth = 0.5, color = "grey40") +
  scale_fill_manual(values = c("Female" = "hotpink3", "Male" = "steelblue")) +
  geom_vline(
    xintercept = mt_cutoff,
    color      = "dodgerblue",
    linetype   = "dashed",
    linewidth  = 0.8
  ) +
  annotate(
    geom  = "text",
    x     = mt_cutoff + 0.5,
    y     = Inf,
    label = paste0("cutoff = ", mt_cutoff, "%"),
    color = "black",
    hjust = 0,
    vjust = 2,
    size  = 3
  ) +
  facet_wrap(~ sample_id, ncol = 4) +
  labs(
    x     = "% Mitochondrial Reads",
    y     = "# of Cells",
    title = "Mitochondrial Read Distribution per Sample (before filtering)",
    fill  = "Sex"
  ) +
  theme_minimal()

## Reshape filter summary to long format for stacked bar plot
filter_summary_long <- reshape(
  filter_summary,
  varying   = c("Kept", "Removed"),
  v.names   = "n_cells",
  timevar   = "status",
  times     = c("Kept", "Removed"),
  direction = "long"
)

## Stacked bar plot of cells retained vs removed per sample
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
    size  = 3,
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
  labs(
    x     = "Sample",
    y     = "# of Cells",
    title = "Cells Retained vs. Removed After QC Filtering",
    fill  = NULL
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


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

message("QC plots saved to plots/")


# Session info ------------------------------------------------------------

sessionInfo()