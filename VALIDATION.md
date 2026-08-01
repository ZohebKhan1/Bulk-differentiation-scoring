# Validation record

Validation was run on 2026-08-01 with native arm64 R 4.6.1 on macOS 26.5.1. The pre-refactor user state is preserved in commit `ee5b8c4` (`save current tutorial revisions`), and the first macOS refactor is in `24c0a00` (`refactor differentiation scoring tutorial`). External R Codex Utils logs contain the complete command output; this file records the durable evidence needed to interpret the repository.

## Input and interface checks

- The three tracked inputs contain 192 aligned samples; counts and VST contain 13,615 identically ordered genes.
- The full `~ cell_line + day_factor` design has rank 27 of 27.
- All seven maintained/test R files passed isolated parse checks after the second-pass refactor.
- `tests/test_score_differentiation_timing.R` passed sample-alignment, reference/held-out, bounded-score, finite-time, finite-distance, and duplicate-ID contract checks.

Input SHA-256 values are recorded in `DATA.md`.

## Cell-line cross-validation

The canonical leave-one-cell-line-out command completed successfully in 1,063 seconds and produced `tmp/GSE122380_leave_one_cell_line_out_validation.rds`:

- 13 held-out cell lines;
- 39 fitted gene-set models;
- 576 score rows across all temporal, maturation, and progenitor sets; and
- a 1,197,455-byte cache with atomic fold-level checkpoints.

For the primary all-temporal model:

| Metric | Value |
| --- | ---: |
| Pearson correlation | 0.9706872292 |
| Squared Pearson correlation | 0.9422336970 |
| Mean absolute error | 0.7021358814 days |
| Median absolute error | 0.4550976693 days |
| Within 1 day | 79.6875% |
| Within 2 days | 94.7917% |

Comparison with the preceding 1,198,850-byte cache found identical score keys and observed days, identical fold summaries, and identical membership for every temporal and direction-split gene set. The maximum absolute predicted-day and residual differences were both zero. The superseded cache contained only legacy names and redundant constant fields and was removed after this comparison.

## Canonical tutorial render

The final render succeeded in 133 seconds. It rebuilt the tutorial objects in a base-parented isolated environment, rendered the single maintained R Markdown source without changing the working directory, and recreated all ten PNG/SVG figure pairs plus the deployed site. The preceding cold render also refreshed the GO cache successfully.

Against the preceding `24c0a00` site, eight of ten final PNG/SVG pairs are byte-identical. All PNG dimensions and colorspaces match. Side-by-side inspection found the other two pairs materially equivalent in scientific content and design:

- The heatmap uses an explicit Helvetica-compatible metric registration for ComplexHeatmap's internal PDF device. This removes repeated macOS font warnings and causes only a negligible layout shift; the matrix, clusters, colors, dimensions, and typography are unchanged.
- The held-out panel filename and titles now say `cell line` rather than the ambiguous `line`. Its score points and weighted day-mean LOESS curves are numerically unchanged to floating-point precision.

The build now emits one audited message for the three of 1,500 temporal profiles that trigger robust-LOESS conditioning warnings; all smoothed values were verified finite, and the former nonfinite-to-zero fallback was removed. The only remaining R warnings are expected GO Ensembl-to-Entrez mapping notices, each affecting 0.05% to 0.33% of the queried IDs. The cold GO-cache rebuild emitted nine notices; the final cached rerender emitted three.

The optimized SVGs use relative paths to the four deployed Nimbus Sans faces; the font-face family matches emitted SVG text. Local browser validation returned HTTP 200 for the page, all ten figures, supporting libraries, font assets, and search index. The page contains 19 inline math spans and eight display equations. The only 404 was the browser's optional `favicon.ico` request.

## Remaining limits

- Reproducibility begins from processed RDS inputs; upstream FASTQ processing and VST construction are unavailable.
- Cell-line validation shares the tracked VST matrix across folds and is not external-cohort validation.
- Exact package versions are recorded in the rendered session information, but the repository does not include an environment lockfile.
- No software or data license was added because license selection requires the repository owner’s explicit choice. Public visibility alone does not grant reuse rights.
