# Analysis and validation decisions

This file records decisions that materially affect the differentiation timing score and tutorial output. Shared numeric values are maintained in `config/analysis.yml`.

## Temporal-gene selection

For the reference samples in each analysis:

1. edgeR TMM factors are estimated from raw counts.
2. A gene is expression-supported when its maximum day-mean normalized CPM is at least 10.
3. DESeq2 fits `~ cell_line + day_factor` and performs an LRT against `~ cell_line`.
4. Genes require an adjusted LRT p-value below `1e-7`.
5. Genes require a day-mean VST range of at least 0.6.

Day is categorical in the LRT. The test therefore detects arbitrary day-dependent expression after adjustment for cell line; it is not a linear-trend test.

## Timing model

- PCA is trained on reference-sample VST values for all selected temporal genes.
- Genes are centered and not variance-scaled.
- The first three PCs define Euclidean projection space.
- Reference samples are averaged within ordered numeric day.
- Consecutive day centroids form finite segments; zero-length segments are omitted.
- Each sample is assigned to the closest finite-segment projection.
- Predicted time is interpolated within that segment and normalized to the first and last reference day.

Segment fractions are clipped to `[0, 1]`. The scorer therefore returns times and scores bounded by the reference interval. Squared distance to the selected projection is returned as an applicability diagnostic but does not change the score.

## Interpretive analyses

The first 1,500 final temporal genes after ordering by LRT adjusted p-value and then raw p-value are shown in the heatmap and split into four clusters. This display subset, its early/late cluster combinations, loading-gene labels, and GO enrichment are interpretive; the primary score uses the full temporal set. Cluster GO tests use the 1,500 displayed genes as their universe. PC1 loading-direction GO tests use the full temporal set as their universe. Ambiguous Ensembl-to-symbol mappings retain the Ensembl identifier.

## Cell-line cross-validation

Each fold holds out one cell line and refits temporal selection, PCA, and the polyline on the other lines. The count-model LRT is fold-specific. The tracked VST object is shared because upstream construction code is unavailable, so the evaluation is internal cell-line cross-validation conditioned on the processed cohort.

All 13 lines contribute to aggregate metrics. Line `19190` is excluded only from the per-line panel grid. Reported `r²` is the square of Pearson correlation and is labeled as such; it is not predictive \(R^2\). The established per-day error figure labels each box with the lower of mean and median absolute error. This is a nonstandard descriptive annotation retained for output continuity; conventional overall mean and median absolute errors are reported separately.

## Cache and output admission

The validation cache is admitted only when its input MD5 values, implementation MD5 values, thresholds, PC count, relevant package versions, and held-out lines match the current run. The tutorial render fails rather than silently omitting or reusing stale validation results. The deployed site is generated only from `tutorial/tutorial.Rmd`; `docs/` contains no maintained R Markdown or build configuration copies.
