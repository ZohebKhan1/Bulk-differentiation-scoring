# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/select_temporal_genes.R
# ----

#' Select temporal genes from counts and VST expression
#'
#' Genes must pass a maximum day-mean TMM CPM threshold, a DESeq2 likelihood-
#' ratio test for categorical day after adjustment for cell line, and a minimum
#' day-mean VST range. Inputs are expected to come from the validated project
#' loader; this function owns the model and filtering contract.
#'
#' @param counts Integer-like raw-count matrix with genes in rows and samples in
#'   columns.
#' @param vst Numeric variance-stabilized expression matrix aligned to `counts`.
#' @param metadata Sample metadata aligned to the matrix columns, including
#'   `sample_id`, `cell_line`, and numeric `day_numeric` columns.
#' @param expression_cpm_cutoff Minimum maximum day-mean TMM CPM.
#' @param lrt_padj_cutoff Maximum adjusted p-value from the DESeq2 LRT.
#' @param vst_dynamic_range_cutoff Minimum range across day-mean VST values.
#'
#' @return A list containing ordered temporal genes, LRT results, day-mean TMM
#'   CPM and VST matrices, and a one-row filtering summary.
select_temporal_genes <- function(
    counts,
    vst,
    metadata,
    expression_cpm_cutoff,
    lrt_padj_cutoff,
    vst_dynamic_range_cutoff) {
  days <- sort(unique(metadata$day_numeric))
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
  lrt_metadata$day_factor <- factor(lrt_metadata$day_numeric, levels = days)
  lrt_metadata$cell_line <- droplevels(factor(lrt_metadata$cell_line))

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
    input_genes = nrow(counts),
    genes_passing_day_mean_tmm_cpm = length(expression_genes),
    genes_tested_by_lrt = nrow(lrt_results),
    genes_passing_lrt_padj = length(lrt_genes),
    genes_passing_all_filters = length(temporal_genes),
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
