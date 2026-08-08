# Differentiation timing score

Bulk RNA-seq differentiation studies are usually sampled at discrete timepoints, but individual samples can progress at different rates. `score_differentiation_timing()` places each sample along a reference differentiation time course using its expression profile rather than its collection day alone.

The function learns a centered, unscaled PCA space from the mean temporal-gene expression profile at each reference timepoint. It automatically retains the smallest number of PCs explaining at least 99% of the between-timepoint variance, connects successive reference-timepoint centroids, and projects each sample to the nearest point on that finite polyline. The result includes a predicted reference time, a normalized score from 0 at the earliest reference timepoint to 1 at the latest, and the sample's distance from the fitted trajectory.

[Open the rendered tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/).

## Example data used in the workflow/tutorial

The worked example uses processed data derived from public [NCBI GEO accession GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380), a bulk RNA-seq time course of human induced pluripotent stem cell differentiation into cardiomyocytes. The analysis includes 192 samples from 13 cell lines collected across differentiation days 1 through 15. These samples form the reference cohort used for temporal-gene selection, PCA training, day-centroid trajectory fitting, and internal validation.

The repository bundles the processed inputs needed by the tutorial:

| File | Contents |
| --- | --- |
| `data/GSE122380_metadata.rds` | Metadata for 192 samples |
| `data/GSE122380_counts.rds` | Raw gene-count matrix with 13,615 genes and 192 samples |
| `data/GSE122380_vst.rds` | Variance-stabilized expression matrix with the same genes and samples |

Raw-read acquisition, sequencing and sample-level QC, alignment, quantification, upstream gene filtering, and construction of the VST matrix are not included. The worked example is therefore reproducible from these processed inputs rather than from FASTQ files.

## Select temporal genes

`select_temporal_genes()` takes raw integer counts, VST expression, and matching
sample metadata. It calculates TMM CPM internally and tests categorical time
with a DESeq2 likelihood-ratio test. Adjustment covariates are specified as
metadata column names; the function constructs the full and reduced formulas.

```r
source('functions/select_temporal_genes.R')

temporal_selection <- select_temporal_genes(
  raw_counts = counts,
  vst_expression = vst,
  metadata = metadata,
  sample_id_col = 'sample_id',
  time_col = 'day_numeric',
  adjustment_covariates = 'cell_line'
)
```

The defaults are `adjustment_covariates = NULL`, `expression_cpm_cutoff = 10`,
`lrt_padj_cutoff = 1e-7`, and `vst_dynamic_range_cutoff = 0.6`. Provide one or
more metadata columns when the categorical-time LRT should adjust for variables
such as batch or donor. Numeric adjustment columns remain continuous; character
and logical columns are treated as factors. To learn temporal genes from one
cohort, supply both `reference_group_col` and `reference_group_value`, for
example `condition` and `control`. These arguments select which samples enter
gene selection rather than setting a DESeq2 coefficient baseline. The function
checks that requested columns and values exist, that the matrices and metadata
contain matching identifiers, and that the resulting DESeq2 model matrices are
estimable before fitting them.

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
  reference_samples = sample_metadata$sample_id[
    sample_metadata$condition == 'control'
  ]
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

Render the tutorial directly:

```bash
Rscript scripts/03_render_tutorial_site.R
```

The render uses `cache/GSE122380_leave_one_cell_line_out_validation.rds` when a compatible local cache exists. Otherwise it runs the complete leave-one-cell-line-out validation in memory and generates the figures from scratch without saving the validation object. To create the optional ignored cache explicitly before rendering, run:

```bash
Rscript scripts/02_run_leave_one_cell_line_out_validation.R
```

The reusable validation computation is defined in `src/run_leave_one_cell_line_out_validation.R`; the cache-building script is only a convenience entry point. The render writes the site to `docs/index.html`.

## Analysis details

The GSE122380 example retains genes with a maximum day-mean TMM CPM of at least 10, a DESeq2 likelihood-ratio-test adjusted p-value below `1e-7`, and a day-mean VST range of at least 0.6. PCA is fitted to the unscaled mean VST profile at each training timepoint after gene centering. The smallest number of PCs explaining at least 99% of the training between-timepoint variance is retained automatically. Ordered timepoint centroids define the finite scoring polyline.

Leave-one-cell-line-out validation repeats temporal-gene selection, PCA training, and trajectory fitting without the held-out cell line. The tracked VST matrix is shared across folds because its upstream construction pipeline is unavailable. The reported results therefore measure internal performance within this processed cohort, not transferability to another dataset or platform.

The score measures timing relative to the chosen reference trajectory. It does not by itself measure cell quality, differentiation efficiency, potency, or lineage identity. Full methods, validation results, and session information are included in the rendered tutorial.

## Repository contents

```text
├── data/         Processed GSE122380 inputs
├── docs/         Generated GitHub Pages site
├── functions/    Reusable temporal-gene selection and scoring functions
├── scripts/      Analysis, cross-validation, and site rendering
├── src/          Internal data-loading and validation helpers
├── tests/        Focused public-function regression tests
└── tutorial/     Maintained tutorial source and presentation assets
```

`tutorial/tutorial.Rmd` is the only maintained site source; `docs/` contains generated output.

## References

1. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550.
2. Strober BJ, Elorbany R, Rhodes K, et al. Dynamic genetic regulation of gene expression during cellular differentiation. *Science*. 2019;364:1287-1290.
