# 13LGS_PilotAnalyses

scRNA-seq analysis of 13-lined ground squirrel (*Ictidomys tridecemlineatus*) ileum.
8 samples (4 female, 4 male) processed through a two-pass pipeline: first-pass clustering
for cell type annotation, followed by ambient RNA denoising with CellSweep.

## Samples

| Sample ID | Sex |
|---|---|
| 13LGS_sIBA_F1_1 | Female |
| 13LGS_sIBA_F1_2 | Female |
| 13LGS_sIBA_F1_3 | Female |
| 13LGS_sIBA_F1_4 | Female |
| 13LGS_sIBA_M1_1 | Male |
| 13LGS_sIBA_M1_2 | Male |
| 13LGS_sIBA_M1_3 | Male |
| 13LGS_sIBA_M1_4 | Male |

## Pipeline

### First Pass — Clustering and Annotation

1. **`scripts/processing/qc_filtering.R`** — QC filtering per sample (MT%, nFeature, nCount cutoffs) and doublet scoring with [`scDblFinder`](https://github.com/plger/scDblFinder)
2. **`scripts/processing/doublet_removal.R`** — Per-sample doublet removal (SLURM array)
3. **`scripts/processing/sctransform_per_sample.R`** — Per-sample SCTransform normalization (SLURM array; `input_dir: "doublets"`, `output_dir: "sct"`)
4. **`scripts/processing/normalize_cluster.R`** — Merge, PCA, Harmony batch correction, UMAP, and clustering (`input_dir: "sct"`, `output_suffix: ""`)
5. **`scripts/analysis/marker_genes.R`** — `FindAllMarkers` per cluster (SLURM job; `input_suffix: ""`)
6. **`scripts/processing/extract_lightweight_data.R`** — Extract UMAP embeddings, metadata, and marker expression from full Seurat object (SLURM job; `input_suffix: ""`)

### `cellsweep` Denoising

7. **`scripts/processing/cellsweep_prep.R`** — Export per-sample raw counts and first-pass cell type labels for `cellsweep`
8. **`scripts/processing/cellsweep_convert.py`** — Convert per-sample MTX files to `AnnData` `.h5ad` (SLURM array)
9. **`scripts/processing/cellsweep_run.py`** — Run `cellsweep` ambient RNA denoising (SLURM array)
10. **`scripts/processing/cellsweep_reload.R`** — Inject denoised counts back into per-sample `Seurat` objects (SLURM job)

### Second Pass — Denoised Clustering and Visualization

11. **`scripts/processing/sctransform_per_sample.R`** — Re-run SCTransform on denoised counts (SLURM array; `input_dir: "cellsweep_clean"`, `output_dir: "sct_denoised"`)
12. **`scripts/processing/normalize_cluster.R`** — Re-run PCA, Harmony, UMAP, and clustering on denoised data (`input_dir: "sct_denoised"`, `output_suffix: "_denoised"`)
13. **`scripts/analysis/marker_genes.R`** — `FindAllMarkers` on denoised clusters (SLURM job; `input_suffix: "_denoised"`)
14. **`scripts/processing/extract_lightweight_data.R`** — Extract lightweight files from denoised `Seurat` object (SLURM job; `input_suffix: "_denoised"`)

### Visualization

15. **`scripts/plots/umap_plots.R`** — UMAP plots colored by sample, sex, cluster, MT%, and cell type (`input_suffix: ""` and `"_denoised"`)
16. **`scripts/plots/marker_plots.R`** — Dot plot, heatmap, and volcano plots for marker genes (`input_suffix: ""` and `"_denoised"`)
17. **`scripts/plots/qc_plots.R`** — QC violin plots and filter summary

## Directory Structure

    data/
      raw/                            — CellRanger .h5 output (filtered + raw)
      processed/
        doublets/                     — Per-sample doublet-cleaned Seurat objects
        sct/                          — Per-sample SCTransform objects (first pass)
        seurat_processed.rds          — Merged, clustered Seurat object (first pass)
        cellsweep/                    — CellSweep input/output per sample
        cellsweep_clean/              — Per-sample Seurat objects with denoised counts
        sct_denoised/                 — Per-sample SCTransform objects (denoised)
        seurat_processed_denoised.rds — Merged, clustered Seurat object (denoised)
    plots/                            — Output figures (.pdf, gitignored)
    scripts/
      processing/                     — Data processing scripts
      analysis/                       — Analysis scripts
      plots/                          — Visualization scripts
      slurm/                          — SLURM submission scripts

## Dependencies

R packages: `Seurat`, `harmony`, `scDblFinder`, `zellkonverter`, `SingleCellExperiment`,
`Matrix`, `ggplot2`, `ggrepel`, `patchwork`, `here`

Python packages: `cellsweep`, `anndata`, `scanpy`, `scipy`, `pandas`, `numpy`

Package versions are managed with `renv` (see `renv.lock`).
