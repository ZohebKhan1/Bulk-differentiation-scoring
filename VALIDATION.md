# Validation record

Validation was run on 2026-08-01 with native arm64 R 4.6.1 on macOS 26.5.1. The pre-refactor user state is preserved in commit `ee5b8c4` (`save current tutorial revisions`). External R Codex Utils logs contain the complete command output; this file records the durable evidence needed to interpret the repository.

## Input and interface checks

- The three tracked inputs contain 192 aligned samples; counts and VST contain 13,615 identically ordered genes.
- The full `~ cell_line + day_factor` design has rank 27 of 27.
- All seven maintained/test R files passed isolated parse checks.
- `tests/test_score_differentiation_timing.R` passed sample-alignment, reference/held-out, bounded-score, finite-distance, and duplicate-ID contract checks.

Input SHA-256 values are recorded in `DATA.md`.

## Cell-line cross-validation

The canonical validation command completed successfully in 1,046 seconds and produced `tmp/GSE122380_leave_one_line_out_validation.rds`:

- 13 held-out cell lines;
- 39 fitted gene-set models;
- 576 score rows across all temporal, maturation, and progenitor sets; and
- a 1,198,850-byte cache, reduced from approximately 17 MB across the two obsolete legacy caches.

For the primary all-temporal model:

| Metric | Value |
| --- | ---: |
| Pearson correlation | 0.9706872292 |
| Squared Pearson correlation | 0.9422336970 |
| Mean absolute error | 0.7021358814 days |
| Median absolute error | 0.4550976693 days |
| Within 1 day | 79.6875% |
| Within 2 days | 94.7917% |

Comparison with the legacy trajectory cache found identical score keys, identical summary tables, and identical membership for every temporal and direction-split gene set. The maximum absolute predicted-day and residual differences were both `3.1974423109204508e-14`, attributable to harmless feature-order floating-point effects.

## Canonical tutorial render

The complete render succeeded in 145 seconds after resolving an initial macOS system-font registration collision. It recreated all 24 admitted outputs: `docs/index.html`, deployed CSS and marker files, the GO cache, and ten PNG/SVG figure pairs.

All PNG dimensions and colorspaces match the baseline. Four PNGs are byte-identical to `ee5b8c4`: the dataset overview, temporal-cluster panel, timing polyline, and per-day LOO error panel. Side-by-side inspection found the other figures materially equivalent in data, geometry, and visual design. Intentional visible differences are limited to:

- `r²` and “Squared Pearson correlation” labels replacing ambiguous `R²`/“R-squared” wording;
- “Differentiation timing score” replacing “Maturation score” on the score figure; and
- retention of an Ensembl ID for one ambiguous symbol mapping.

The optimized SVGs use relative paths to the four deployed Nimbus Sans faces; the font-face family matches emitted SVG text. Local browser validation returned HTTP 200, loaded all ten figures without broken assets, found 27 math spans and eight display equations, activated three code-fold sections, loaded the bundled Latin Modern Sans site face, and retrieved the generated search index. The only 404 was the browser’s optional `favicon.ico` request.

## Remaining limits

- Reproducibility begins from processed RDS inputs; upstream FASTQ processing and VST construction are unavailable.
- Cell-line validation shares the tracked VST matrix across folds and is not external-cohort validation.
- Exact package versions are recorded in the rendered session information, but the repository does not include an environment lockfile.
- No software or data license was added because license selection requires the repository owner’s explicit choice. Public visibility alone does not grant reuse rights.
