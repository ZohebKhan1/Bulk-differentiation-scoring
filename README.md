# Differentiation timing score

Samples collected on the same nominal differentiation day can be at different molecular stages. A collection date therefore provides an experimental label, but not always a comparable measure of progress.

This workflow uses an ordered reference time course to estimate each sample's position along a molecular differentiation trajectory. Temporal genes define a centered, unscaled PCA space; ordered reference-timepoint centroids define a finite trajectory through that space; and each sample is projected to its nearest point on the trajectory. The result is a normalized predicted reference time, a differentiation score bounded by the reference endpoints, and a squared projection distance that indicates how well the sample is represented by that trajectory.

## When to use this workflow

Use the workflow when you have a suitable ordered reference differentiation time course and want to compare molecular timing within that reference frame. The analysis requires:

- normalized genes-by-samples expression with stable gene and sample identifiers;
- sample metadata aligned to the expression columns, including ordered reference times; and
- temporal genes selected independently of the samples that will be scored.

The score is relative to the chosen reference and its normalization. It is not a cell-quality metric, potency or differentiation-efficiency measure, lineage-identity classifier, or estimate of absolute biological age. Samples far from the fitted trajectory should be treated as model-relative diagnostics rather than forced into a biological interpretation.

## Start here

[Read the rendered tutorial](https://zohebkhan1.github.io/pca-maturation-scoring/) for the canonical explanation, figures, validation results, worked example, and interpretation guidance. The repository also contains reusable implementations for temporal-gene selection and differentiation-timing scoring in `functions/`.

The returned predicted reference time stays on the original reference-time scale, while the normalized score makes endpoint-scaled comparisons within that same reference easier to read. The squared distance is useful for finding samples that do not resemble the learned trajectory. These outputs should be reviewed together: a timing estimate without its reference context is not an absolute measurement.

## Case study and reproducibility scope

The tutorial uses processed inputs derived from [NCBI GEO accession GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380), a bulk RNA-seq time course of human induced pluripotent stem cell differentiation into cardiomyocytes. Metadata, raw counts, and variance-stabilized expression are provided so the documented workflow can be reproduced from processed data. Raw FASTQ acquisition, alignment, quantification, upstream filtering, and construction of the VST matrix are outside this repository's scope.

Leave-one-cell-line-out validation measures performance within this processed cohort. It does not establish that the score will transfer unchanged to another study, platform, protocol, or reference time course; adaptation requires a reference appropriate to the new setting.

## Repository layout

```text
├── data/         Processed GSE122380 inputs
├── docs/         Rendered GitHub Pages site
├── functions/    Reusable temporal-gene selection and scoring functions
├── scripts/      Analysis and site-rendering scripts
├── src/          Internal data-loading and validation helpers
├── tests/        Public-function regression tests
└── tutorial/     Maintained tutorial source and presentation assets
```

## References

1. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550.
2. Strober BJ, Elorbany R, Rhodes K, et al. Dynamic genetic regulation of gene expression during cellular differentiation. *Science*. 2019;364:1287-1290.
