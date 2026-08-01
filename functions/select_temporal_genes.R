# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/select_temporal_genes.R
# ----

# Select temporal genes with day-level expression, LRT, and VST-range evidence
select_temporal_genes <- function(
    counts,
    vst,
    metadata,
    expression_cpm_cutoff,
    lrt_padj_cutoff,
    vst_dynamic_range_cutoff) {
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop('counts must be a numeric matrix.', call. = FALSE)
  }
  if (!is.matrix(vst) || !is.numeric(vst)) {
    stop('vst must be a numeric matrix.', call. = FALSE)
  }
  if (
    is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))
  ) {
    stop('counts must have unique gene and sample names.', call. = FALSE)
  }
  if (
    !identical(rownames(counts), rownames(vst)) ||
      !identical(colnames(counts), colnames(vst))
  ) {
    stop('counts and vst must have identical ordered genes and samples.', call. = FALSE)
  }
  if (!is.data.frame(metadata)) {
    stop('metadata must be a data frame.', call. = FALSE)
  }
  required_metadata_columns <- c('sample_id', 'day_numeric', 'cell_line')
  missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
  if (length(missing_metadata_columns) > 0L) {
    stop(
      'metadata is missing columns: ',
      paste(missing_metadata_columns, collapse = ', '),
      '.',
      call. = FALSE
    )
  }
  metadata_sample_ids <- as.character(metadata$sample_id)
  if (anyNA(metadata_sample_ids) || anyDuplicated(metadata_sample_ids)) {
    stop('metadata sample IDs must be complete and unique.', call. = FALSE)
  }
  metadata_rows <- match(colnames(counts), metadata_sample_ids)
  if (anyNA(metadata_rows) || nrow(metadata) != ncol(counts)) {
    stop('metadata must contain exactly one row for every matrix sample.', call. = FALSE)
  }
  metadata <- metadata[metadata_rows, , drop = FALSE]
  if (!is.numeric(metadata$day_numeric) || any(!is.finite(metadata$day_numeric))) {
    stop('metadata day_numeric values must be finite numeric values.', call. = FALSE)
  }

  thresholds <- c(
    expression_cpm_cutoff = expression_cpm_cutoff,
    lrt_padj_cutoff = lrt_padj_cutoff,
    vst_dynamic_range_cutoff = vst_dynamic_range_cutoff
  )
  if (
    !is.numeric(thresholds) || anyNA(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0)
  ) {
    stop('Temporal-gene thresholds must be finite nonnegative numbers.', call. = FALSE)
  }

  days <- sort(unique(metadata$day_numeric))
  if (length(days) < 2L) {
    stop('metadata must span at least two differentiation days.', call. = FALSE)
  }
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::normLibSizes(dge, method = 'TMM')
  tmm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  day_mean_tmm_cpm <- sapply(days, function(day_value) {
    day_samples <- metadata$sample_id[metadata$day_numeric == day_value]
    rowMeans(tmm_cpm[, day_samples, drop = FALSE], na.rm = TRUE)
  })
  colnames(day_mean_tmm_cpm) <- paste0('D', days)
  max_day_mean_tmm_cpm <- apply(day_mean_tmm_cpm, 1, max, na.rm = TRUE)
  expression_genes <- names(max_day_mean_tmm_cpm)[
    max_day_mean_tmm_cpm >= expression_cpm_cutoff
  ]
  if (length(expression_genes) == 0L) {
    stop('No genes passed the day-mean TMM CPM threshold.', call. = FALSE)
  }

  lrt_metadata <- metadata
  rownames(lrt_metadata) <- lrt_metadata$sample_id
  lrt_metadata$day_factor <- droplevels(lrt_metadata$day_factor)
  lrt_metadata$cell_line <- droplevels(lrt_metadata$cell_line)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts[expression_genes, , drop = FALSE],
    colData = lrt_metadata,
    design = ~ cell_line + day_factor
  )
  dds <- DESeq2::DESeq(
    dds,
    test = 'LRT',
    reduced = ~ cell_line,
    quiet = TRUE
  )

  lrt_results <- as.data.frame(DESeq2::results(dds, alpha = lrt_padj_cutoff))
  lrt_results$gene_id <- rownames(lrt_results)
  lrt_results <- lrt_results[
    order(lrt_results$padj, lrt_results$pvalue, na.last = TRUE),
    ,
    drop = FALSE
  ]
  lrt_genes <- lrt_results$gene_id[
    !is.na(lrt_results$padj) & lrt_results$padj < lrt_padj_cutoff
  ]
  if (length(lrt_genes) == 0L) {
    stop('No genes passed the DESeq2 LRT adjusted p-value threshold.', call. = FALSE)
  }

  day_mean_vst <- sapply(days, function(day_value) {
    day_samples <- metadata$sample_id[metadata$day_numeric == day_value]
    rowMeans(vst[lrt_genes, day_samples, drop = FALSE], na.rm = TRUE)
  })
  colnames(day_mean_vst) <- paste0('D', days)
  vst_dynamic_range <- apply(day_mean_vst, 1, function(gene_values) {
    max(gene_values, na.rm = TRUE) - min(gene_values, na.rm = TRUE)
  })
  temporal_genes <- lrt_genes[
    vst_dynamic_range[lrt_genes] >= vst_dynamic_range_cutoff
  ]
  if (length(temporal_genes) < 2L) {
    stop('Fewer than two genes passed all temporal-gene filters.', call. = FALSE)
  }

  summary <- data.frame(
    total_expressed_genes = nrow(counts),
    genes_passing_day_mean_tmm_cpm = length(expression_genes),
    genes_tested_by_lrt = nrow(lrt_results),
    genes_passing_lrt_padj = length(lrt_genes),
    genes_passing_vst_dynamic_range = length(temporal_genes),
    expression_cpm_cutoff = expression_cpm_cutoff,
    lrt_padj_cutoff = lrt_padj_cutoff,
    vst_dynamic_range_cutoff = vst_dynamic_range_cutoff,
    stringsAsFactors = FALSE
  )

  list(
    temporal_genes = temporal_genes,
    lrt_results = lrt_results,
    day_mean_tmm_cpm = day_mean_tmm_cpm,
    day_mean_vst = day_mean_vst,
    summary = summary
  )
}
