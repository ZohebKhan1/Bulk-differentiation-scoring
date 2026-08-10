# Differentiation timing score

This repository provides R functions that select genes whose expression changes across an ordered reference time course and project samples onto the resulting differentiation trajectory. Given raw counts, variance-stabilized expression, and aligned metadata, the workflow returns a relative predicted reference time, an endpoint-scaled differentiation score, and a squared projection distance.

The [rendered tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/) is the canonical guide to the method, worked example, figures, validation, and adaptation. This README focuses on the functions and the smallest reproducible start.

## Functions

| Function | Use | Main result |
| --- | --- | --- |
| `select_temporal_genes()` | Filter a reference cohort using maximum day-mean TMM CPM, a categorical-time DESeq2 LRT, and minimum day-mean VST range. | Ordered temporal-gene IDs, LRT results, and a filtering summary. |
| `score_differentiation_timing()` | Fit a reference trajectory and project samples onto it. | Sample-level timing scores, PCA coordinates, the centroid polyline, and projection diagnostics. |

The scorer fits centered, unscaled PCA to the mean expression profile at each reference timepoint. It connects ordered timepoint centroids with finite line segments and projects each sample to its nearest segment. The function retains enough PCs to explain at least 99% of the between-timepoint variance.

## Quick start

The repository includes processed inputs from [GEO accession GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380), a bulk RNA-seq time course of human induced pluripotent stem cell differentiation into cardiomyocytes. The example contains 192 samples from 13 cell lines across 15 differentiation days.

Temporal-gene selection requires `DESeq2` and `edgeR` from Bioconductor. Install them before running the example:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
```

Run from the repository root:

```r
source("functions/select_temporal_genes.R")
source("functions/score_differentiation_timing.R")

metadata <- readRDS("data/GSE122380_metadata.rds")
counts <- readRDS("data/GSE122380_counts.rds")
vst <- readRDS("data/GSE122380_vst.rds")

# Keep both matrices in the metadata sample order.
sample_ids <- metadata$sample_id
counts <- counts[, sample_ids, drop = FALSE]
vst <- vst[, sample_ids, drop = FALSE]

temporal_selection <- select_temporal_genes(
  raw_counts = counts,
  vst_expression = vst,
  metadata = metadata,
  adjustment_covariates = "cell_line"
)

timing_fit <- score_differentiation_timing(
  expression_matrix = vst,
  metadata = metadata,
  temporal_genes = temporal_selection$temporal_genes
)

timing_fit$scores[, c(
  "sample_id",
  "observed_time",
  "predicted_time",
  "differentiation_score",
  "squared_distance"
)]
```

This example fits and scores the included reference cohort for calibration. It is not held-out validation. For new samples, select genes and fit the trajectory on independent reference samples. Then supply those IDs through `reference_samples`. Here, `cell_line` adjusts the categorical-time test during gene selection. Choose covariates that match your reference study. The complete workflow and case-study filter settings are documented in the [tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/).

## Input requirements

- `select_temporal_genes()` requires integer-like, unnormalized counts; it calculates TMM CPM internally. It also requires VST expression with the same gene and sample identifiers as the count matrix.
- `select_temporal_genes()` requires sample metadata with one row per expression sample, matching sample IDs, and finite numeric timepoints.
- `score_differentiation_timing()` requires matching sample IDs, normalized expression such as VST expression, and temporal genes selected without using the samples being scored. Reference samples must have finite timepoints; non-reference samples can have missing timepoints. Supply `reference_samples` to fit the PCA and trajectory on a reference cohort before scoring new samples.

## Interpret the result

`timing_fit$scores` reports one row per sample:

- `predicted_time` is the estimated position on the original reference-time scale.
- `differentiation_score` maps the earliest and latest reference timepoints to 0 and 1, respectively.
- `squared_distance` measures how far the sample lies from its nearest fitted trajectory segment in PCA space.

The score is relative to the selected reference. Compare scores within the same reference frame and inspect `squared_distance` alongside the timing estimate when a sample may not fit the learned trajectory.

## Data and further documentation

The tracked `data/` files are processed starting points for the case study: metadata, raw counts, and VST expression. The repository starts from these processed matrices; see the tutorial for the upstream-processing boundary and assumptions.

Use the [tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/) for the full workflow, figures, leave-one-cell-line-out validation, interpretation, and adaptation guidance. The reusable implementations are [`functions/select_temporal_genes.R`](functions/select_temporal_genes.R) and [`functions/score_differentiation_timing.R`](functions/score_differentiation_timing.R).
