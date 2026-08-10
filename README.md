# Differentiation timing score

[See differentiation maturation scoring example here.](https://zohebkhan1.github.io/pca-maturation-scoring/)

Use this workflow with bulk RNA-seq samples from a differentiation time course to estimate how far each sample has progressed through differentiation. The scoring function defines a reference-specific transcriptomic clock in principal-component space.

This repository provides R functions that select genes whose expression changes across an ordered reference time course and project samples onto the differentiation trajectory. Given raw counts, variance-stabilized expression, and aligned metadata, the workflow returns a predicted time on the reference scale, an endpoint-scaled differentiation score, and a squared projection distance.

## API

- [`select_temporal_genes()`](functions/select_temporal_genes.R) selects genes with time-dependent expression from raw counts, VST expression, and aligned sample metadata.
- [`score_differentiation_timing()`](functions/score_differentiation_timing.R) projects samples onto an ordered reference trajectory and returns predicted time, a normalized score, and projection distance.
- [`run_leave_one_cell_line_out_validation()`](src/run_leave_one_cell_line_out_validation.R) refits temporal-gene selection and trajectory scoring for each held-out cell line and returns fold-level predictions and summaries.

## Quick start

Temporal-gene selection requires `DESeq2` and `edgeR` from Bioconductor. Install them before using the function:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
```

Run the following code in an R session that contains your own `counts`, `vst`, and `metadata` objects. The objects must satisfy the input requirements below.

```r
base_url <- "https://raw.githubusercontent.com/ZohebKhan1/pca-maturation-scoring/main/functions"

download.file(
  paste0(base_url, "/select_temporal_genes.R"),
  "select_temporal_genes.R"
)

download.file(
  paste0(base_url, "/score_differentiation_timing.R"),
  "score_differentiation_timing.R"
)

source("select_temporal_genes.R")
source("score_differentiation_timing.R")
```

## Select temporal/time-dependent genes

`select_temporal_genes()` selects genes that change across the reference time course. It applies three filters in sequence: an expression filter, a DESeq2 likelihood-ratio test (LRT) for categorical time, and a VST-range filter. Only genes that pass all three filters enter the trajectory model.

```mermaid
flowchart LR
    A[Reference counts and metadata] --> B[1. Expression filter]
    B --> C[2. DESeq2 LRT for time]
    C --> D[3. VST-range filter]
    D --> E[Temporal genes]
```

```r
temporal_selection <- select_temporal_genes(
  raw_counts = counts,
  vst_expression = vst,
  metadata = metadata
)

temporal_genes <- temporal_selection$temporal_genes
```

## Define the reference differentiation trajectory

`score_differentiation_timing()` fits centered, unscaled PCA to the mean temporal-gene profile at each reference timepoint. It connects the ordered timepoint centroids with finite line segments and projects each sample to its nearest point on the trajectory.

```r
timing_fit <- score_differentiation_timing(
  expression_matrix = vst,
  metadata = metadata,
  temporal_genes = temporal_genes
)

timing_fit$scores[, c(
  "sample_id",
  "observed_time",
  "predicted_time",
  "differentiation_score",
  "squared_distance"
)]
```

The result reports predicted time on the reference scale, an endpoint-scaled score, and distance from the fitted trajectory. For new samples, select genes and fit the trajectory using independent reference samples, then supply their IDs through `reference_samples`.

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
