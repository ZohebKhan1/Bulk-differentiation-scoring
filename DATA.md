# Data contract

The repository tracks the processed inputs required by the tutorial. All three objects must contain the same 192 uniquely identified samples in metadata order. The count and VST matrices must also contain the same 13,615 uniquely identified genes in identical order.

| File | Object contract | Size | SHA-256 |
| --- | --- | ---: | --- |
| `data/GSE122380_metadata.rds` | data frame; 192 samples | 4,944 bytes | `9d65f5d33cdbc9755a3e2222cf24b3a06f263954ed428dfd836012d0679f4279` |
| `data/GSE122380_counts.rds` | integer-valued matrix; 13,615 × 192 | 4,119,440 bytes | `373a1f784477337aa5ec5c21f367e0997719feccaae92c94d490fd8a1cb1cff1` |
| `data/GSE122380_vst.rds` | finite numeric matrix; 13,615 × 192 | 17,656,288 bytes | `fe14d5f200aecc0bf57c41a78ab8c7c0f7de544a75c1bb915cb6d4f0c18951bc` |

`functions/load_GSE122380_data.R` enforces the object types, identifiers, ordering, completeness, finiteness, and nonnegative integer count contract before analysis.

## Provenance boundary

The samples originate from the GSE122380 iPSC-to-cardiomyocyte bulk RNA-seq time course. The repository does not contain the upstream raw-read acquisition, QC, alignment, quantification, gene-filtering, or VST-construction pipeline. Consequently:

- the versioned RDS files are the reproducible analysis boundary;
- the exact checksums above identify the inputs used for the rendered result;
- claims about upstream processing cannot be independently reconstructed here; and
- the shared VST object is used in every cell-line cross-validation fold.

Raw and intermediate data locations are ignored deliberately; tracked processed inputs are not regenerated or overwritten by any maintained script.
