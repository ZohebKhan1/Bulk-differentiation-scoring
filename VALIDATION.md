# Validation record

Validation was run on 2026-08-05 with native arm64 R 4.6.1 on macOS
26.5.1. Input SHA-256 values are recorded in `DATA.md`, and exact package
versions appear in the rendered tutorial's session information.

The maintained checks can be run from the repository root:

```sh
Rscript scripts/02_run_leave_one_cell_line_out_validation.R
Rscript scripts/03_render_tutorial_site.R
```

## Input and interface checks

- The three tracked inputs contain 192 aligned samples; counts and VST contain 13,615 identically ordered genes.
- The full `~ cell_line + day_factor` design has rank 27 of 27.
- All six maintained R files passed isolated parse checks.

## Cell-line cross-validation

The leave-one-cell-line-out validation completed successfully in 1,572 seconds:

- 13 held-out cell lines;
- 39 fitted gene-set models;
- 576 score rows across all temporal, maturation, and progenitor sets.

For the primary all-temporal model:

| Metric | Value |
| --- | ---: |
| Pearson correlation | 0.9706872292 |
| Squared Pearson correlation | 0.9422336970 |
| Mean absolute error | 0.7021358814 days |
| Median absolute error | 0.4550976693 days |
| Within 1 day | 79.6875% |
| Within 2 days | 94.7917% |

Each fold refits temporal-gene selection, PCA, and the scoring polyline without
using the held-out cell line during model fitting.

## Tutorial render

The tutorial render completed in 235 seconds and produced ten high-resolution
PNG figures plus the deployed site. The rendered figures retained the expected
data, dimensions, colors, and typography. Site assets, fonts, local links, and
the search index resolved correctly.

Three of 1,500 temporal profiles triggered robust-LOESS conditioning warnings;
all resulting smoothed values were finite. Expected Ensembl-to-Entrez mapping
notices affected 0.05% to 0.33% of identifiers in the corresponding queries.

## Remaining limits

- Reproducibility begins from processed RDS inputs; upstream FASTQ processing and VST construction are unavailable.
- Cell-line validation shares the tracked VST matrix across folds and is not external-cohort validation.
- Exact package versions are recorded in the rendered session information, but the repository does not include an environment lockfile.
- No software or data license was added because license selection requires the repository owner’s explicit choice. Public visibility alone does not grant reuse rights.
