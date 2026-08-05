# Differentiation timing score

Bulk RNA-seq differentiation studies are usually sampled at discrete timepoints, but individual samples can progress at different rates. `score_differentiation_timing()` places each sample along a reference differentiation time course using its expression profile rather than its collection day alone.

The function learns a centered, unscaled PCA space from reference samples and preselected temporal genes. It connects the average position of successive reference timepoints, then projects each sample to the nearest point on that trajectory. The result includes a predicted reference time, a normalized score from 0 at the earliest reference timepoint to 1 at the latest, and the sample's distance from the fitted trajectory.

The worked example uses 192 processed GSE122380 iPSC-to-cardiomyocyte samples from 13 cell lines across days 1 through 15. [Open the rendered tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/).

## Use the scorer

The scorer is a single dependency-free R file. It expects normalized expression, matching sample metadata, and a temporal-gene set selected independently of the samples being evaluated. It does not normalize counts, select genes, draw figures, or write files.

```r
base_url <- paste0(
  'https://raw.githubusercontent.com/ZohebKhan1/',
  'pca-maturation-scoring/main/functions'
)

download.file(
  paste0(base_url, '/score_differentiation_timing.R'),
  'score_differentiation_timing.R'
)
source('score_differentiation_timing.R')

timing_fit <- score_differentiation_timing(
  expression_matrix = vst_matrix,
  metadata = sample_metadata,
  temporal_genes = temporal_genes,
  sample_id_col = 'sample_id',
  time_col = 'day_numeric',
  reference_col = 'condition',
  reference_values = 'control',
  n_pcs = 3
)

head(timing_fit$scores)
```

`timing_fit$scores` contains the observed and predicted time, the normalized `differentiation_score`, the nearest trajectory segment, and the squared projection distance for each sample. A large distance indicates that a sample is poorly represented by the reference trajectory and should be interpreted cautiously.

## Reproduce the worked example

Use a current R installation. On Apple silicon, install packages for the native arm64 R architecture.

```r
install.packages(c(
  'BiocManager', 'bookdown', 'circlize', 'ggplot2', 'ggrepel',
  'patchwork', 'ragg', 'scales', 'systemfonts', 'viridis'
))

BiocManager::install(c(
  'AnnotationDbi', 'clusterProfiler', 'ComplexHeatmap', 'DESeq2',
  'edgeR', 'org.Hs.eg.db'
))
```

Run the held-out validation before rendering the tutorial:

```bash
Rscript scripts/02_run_leave_one_cell_line_out_validation.R
Rscript scripts/03_render_tutorial_site.R
```

The validation script checkpoints completed folds in `tmp/`. The render rejects a missing or stale cache, rebuilds the analysis figures, and writes the site to `docs/index.html`.

## Analysis details

The GSE122380 example retains genes with a maximum day-mean TMM CPM of at least 10, a DESeq2 likelihood-ratio-test adjusted p-value below `1e-7`, and a day-mean VST range of at least 0.6. PCA is centered without variance scaling and retains three PCs. Ordered day centroids define the finite scoring trajectory.

Leave-one-cell-line-out validation repeats temporal-gene selection, PCA training, and trajectory fitting without the held-out cell line. The tracked VST matrix is shared across folds because its upstream construction pipeline is unavailable. The reported results therefore measure internal performance within this processed cohort, not transferability to another dataset or platform.

See [METHODS.md](METHODS.md) for the analysis decisions, [DATA.md](DATA.md) for input dimensions and checksums, and [VALIDATION.md](VALIDATION.md) for executed checks and remaining limitations.

## Repository contents

```text
├── data/         Processed GSE122380 inputs
├── docs/         Generated GitHub Pages site
├── functions/    Data loading, temporal selection, and scoring functions
├── scripts/      Analysis, cross-validation, and site rendering
└── tutorial/     Maintained tutorial source and presentation assets
```

`tutorial/tutorial.Rmd` is the only maintained site source; `docs/` contains generated output.

## Scope

Reproducibility begins with the three versioned processed RDS objects in `data/`. The repository does not include raw-read acquisition, QC, alignment, quantification, or VST-construction code. The score measures timing relative to the chosen reference trajectory; it does not by itself measure cell quality, differentiation efficiency, potency, or lineage identity.

## References

1. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550.
2. Strober BJ, Elorbany R, Rhodes K, et al. Dynamic genetic regulation of gene expression during cellular differentiation. *Science*. 2019;364:1287-1290.
