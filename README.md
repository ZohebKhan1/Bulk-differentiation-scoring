# Differentiation timing score

[See differentiation maturation scoring example here.](https://zohebkhan1.github.io/pca-maturation-scoring/)

Use this workflow with bulk RNA-seq samples from a differentiation time course to estimate how far each sample has progressed through differentiation. The scoring function defines a reference-specific transcriptomic clock in principal-component space.

This repository provides R functions that select genes whose expression changes across an ordered reference time course and project samples onto the differentiation trajectory. Given raw counts, variance-stabilized expression, and aligned metadata, the workflow returns a predicted time on the reference scale, an endpoint-scaled differentiation score, and a squared projection distance.

## API

- [`select_temporal_genes()`](functions/select_temporal_genes.R) selects genes with time-dependent expression from raw counts, VST expression, and aligned sample metadata.
- [`score_differentiation_timing()`](functions/score_differentiation_timing.R) projects samples onto an ordered reference trajectory and returns predicted time, a normalized score, and projection distance.
- [`run_leave_one_cell_line_out_validation()`](src/run_leave_one_cell_line_out_validation.R) refits temporal-gene selection and trajectory scoring for each held-out cell line and returns fold-level predictions and summaries.

The scorer fits centered, unscaled PCA to the mean expression profile at each reference timepoint. It connects ordered timepoint centroids with finite line segments and projects each sample to its nearest segment. The function retains enough PCs to explain at least 99% of the between-timepoint variance.

## Quick start

Temporal-gene selection requires `DESeq2` and `edgeR` from Bioconductor. Install them before running the example:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
```

Run the following code in an R session that contains your own `counts`, `vst`, and `metadata` objects. The objects must satisfy the input requirements below.

```r
function_base_url <- paste0(
  "https://raw.githubusercontent.com/ZohebKhan1/",
  "pca-maturation-scoring/main/functions"
)

function_files <- c(
  "select_temporal_genes.R",
  "score_differentiation_timing.R"
)

for (function_file in function_files) {
  download.file(
    url = paste(function_base_url, function_file, sep = "/"),
    destfile = function_file,
    mode = "wb"
  )
  source(function_file)
}

temporal_selection <- select_temporal_genes(
  raw_counts = counts,
  vst_expression = vst,
  metadata = metadata
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

The code fits and scores the samples supplied in `metadata` as one reference cohort. For new samples, select genes and fit the trajectory on independent reference samples. Then supply those IDs through `reference_samples`. Add `adjustment_covariates` when your study design includes variables such as batch, donor, or cell line.

## Input requirements

- Both functions use `sample_id` and `day_numeric` as the default metadata columns; set `sample_id_col` and `time_col` when your metadata uses different names.
- `select_temporal_genes()` requires integer-like, unnormalized counts; it calculates TMM CPM internally. It also requires VST expression with the same gene and sample identifiers as the count matrix.
- `select_temporal_genes()` requires sample metadata with one row per expression sample, matching sample IDs, and finite numeric timepoints.
- `score_differentiation_timing()` requires matching sample IDs, normalized expression such as VST expression, and temporal genes selected without using the samples being scored. Reference samples must have finite timepoints; non-reference samples can have missing timepoints. Supply `reference_samples` to fit the PCA and trajectory on a reference cohort before scoring new samples.

## Interpret the result

`timing_fit$scores` reports one row per sample:

- `predicted_time` is the estimated position on the original reference-time scale.
- `differentiation_score` maps the earliest and latest reference timepoints to 0 and 1, respectively.
- `squared_distance` measures how far the sample lies from its nearest fitted trajectory segment in PCA space.

The score is relative to the selected reference. Compare scores within the same reference frame and inspect `squared_distance` alongside the timing estimate when a sample may not fit the learned trajectory.

## Data and tutorial

The repository includes processed metadata, raw counts, and VST expression from a [GEO accession GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380) case study. These files support the tutorial but are not required for the Quick start. The repository starts from processed matrices; see the tutorial for data-preparation details.

Use the [tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/) for the analysis workflow, figures, leave-one-cell-line-out validation, interpretation, and adaptation. The functions are [`functions/select_temporal_genes.R`](functions/select_temporal_genes.R) and [`functions/score_differentiation_timing.R`](functions/score_differentiation_timing.R).

## Citations

1. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
2. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26(1):139–140. [doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)
3. Strober BJ, Elorbany R, Rhodes K, et al. Dynamic genetic regulation of gene expression during cellular differentiation. *Science*. 2019;364(6447):1287–1290. [doi:10.1126/science.aaw0040](https://doi.org/10.1126/science.aaw0040)
4. [NCBI GEO accession GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380)
