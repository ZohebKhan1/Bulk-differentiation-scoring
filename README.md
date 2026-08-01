# Differentiation timing score

A reproducible bulk RNA-seq tutorial and standalone R scorer for estimating position along a reference differentiation time course. The worked example uses 192 processed GSE122380 iPSC-to-cardiomyocyte samples from 13 cell lines across days 1 through 15.

[Open the rendered tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/)

## Reuse the scorer

The scorer is one dependency-free R file. It accepts normalized expression, metadata, and a preselected temporal-gene set; it does not normalize counts, select genes, draw figures, or write files.

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

The function trains PCA only on reference samples, joins ordered reference-time centroids, and projects every sample to its nearest position on the finite centroid polyline. `predicted_time` and `differentiation_score` are bounded by the reference interval. The returned squared distance is useful for identifying samples far from the fitted trajectory.

## Repository layout

```text
├── config/       Shared scientific thresholds and retained-PC count
├── data/         Tracked processed GSE122380 inputs
├── docs/         Generated GitHub Pages site; do not edit by hand
├── functions/    Data boundary, temporal selection, and standalone scorer
├── scripts/      Analysis, cell-line cross-validation, and site render
├── tests/        Dependency-free scorer contract checks
└── tutorial/     Single maintained tutorial source, CSS, and font assets
```

`tutorial/tutorial.Rmd` is the only maintained site source. `scripts/03_render_tutorial_site.R` renders it directly to `docs/`, preventing source and deployed copies from drifting.

## Reproduce on macOS

Use a current R installation. On Apple silicon, all packages should be installed for the native arm64 R architecture. The following separates CRAN and Bioconductor dependencies:

```r
install.packages(c(
  'BiocManager', 'bookdown', 'circlize', 'ggplot2', 'ggrepel',
  'patchwork', 'ragg', 'scales', 'svglite', 'systemfonts',
  'viridis', 'yaml'
))

BiocManager::install(c(
  'AnnotationDbi', 'clusterProfiler', 'ComplexHeatmap', 'DESeq2',
  'edgeR', 'org.Hs.eg.db'
))
```

From the repository root, run the held-out validation first, then render:

```bash
Rscript scripts/02_run_leave_one_cell_line_out_validation.R
Rscript scripts/03_render_tutorial_site.R
```

The validation script checkpoints completed folds in `tmp/`. Its cache key includes input, implementation, parameter, package-version, and held-out-cell-line identity. The render refuses a missing or stale cache, rebuilds the tutorial objects and all ten figure pairs, and writes `docs/index.html`.

Run the fast standalone-scorer checks independently with:

```bash
Rscript tests/test_score_differentiation_timing.R
```

## Scientific contract

Temporal genes must satisfy all three configured criteria: maximum day-mean TMM CPM at least 10, DESeq2 LRT adjusted p-value below `1e-7`, and day-mean VST range at least 0.6. PCA is centered without variance scaling and retains three PCs. Ordered day centroids form the scoring polyline; the day 1 and day 15 centroids anchor scores 0 and 1.

Leave-one-cell-line-out validation refits temporal-gene selection, PCA, and the polyline within each fold. The tracked VST matrix is shared across folds because its upstream construction pipeline is unavailable. These results are internal validation conditional on the processed GSE122380 cohort, not an external transferability estimate.

See [METHODS.md](METHODS.md) for analysis decisions, [DATA.md](DATA.md) for input dimensions and checksums, and [VALIDATION.md](VALIDATION.md) for the executed checks and remaining limits. Exact R and package versions are recorded at the end of the rendered tutorial.

## Scope and provenance

The repository starts from three tracked processed RDS objects. It does not include raw-read acquisition, QC, alignment, quantification, or VST-construction code, so reproducibility begins at these versioned inputs rather than FASTQ files. The score measures trajectory position; it should not be interpreted by itself as cell quality, differentiation efficiency, potency, or lineage identity.

## References

1. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550.
2. Strober BJ, Elorbany R, Rhodes K, et al. Dynamic genetic regulation of gene expression during cellular differentiation. *Science*. 2019;364:1287-1290.
3. Xie Y. *bookdown: Authoring Books and Technical Documents with R Markdown*. 2016.

## Contact

Zoheb Khan, Moskowitz Lab, University of Chicago

zohebkhan600@gmail.com · [zohebkhan1.github.io](https://zohebkhan1.github.io/)
