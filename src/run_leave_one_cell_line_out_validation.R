.leave_one_cell_line_out_validation_settings <- function(
    expression_cpm_cutoff,
    lrt_padj_cutoff,
    vst_dynamic_range_cutoff) {
  list(
    method = 'unscaled_day_mean_pca_99_polyline_v1',
    expression_cpm_cutoff = expression_cpm_cutoff,
    lrt_padj_cutoff = lrt_padj_cutoff,
    vst_dynamic_range_cutoff = vst_dynamic_range_cutoff
  )
}

#' Run leave-one-cell-line-out timing validation
#'
#' Temporal genes, the PCA reference space, and the scoring polyline are refit
#' inside every held-out fold. The function is deliberately free of file I/O so
#' callers can use its return value directly or cache it explicitly.
#'
#' @param counts Raw-count matrix with genes in rows and samples in columns.
#' @param vst Variance-stabilized expression matrix aligned to `counts`.
#' @param metadata Sample metadata aligned to the matrix columns.
#' @param expression_cpm_cutoff Minimum maximum day-mean TMM CPM.
#' @param lrt_padj_cutoff Maximum adjusted p-value from the temporal LRT.
#' @param vst_dynamic_range_cutoff Minimum day-mean VST range.
#'
#' @return A list containing held-out scores, fold summaries, and detailed fold
#'   results. The function does not save the returned object.
run_leave_one_cell_line_out_validation <- function(
    counts,
    vst,
    metadata,
    expression_cpm_cutoff = 10,
    lrt_padj_cutoff = 1e-7,
    vst_dynamic_range_cutoff = 0.6) {
  required_metadata_columns <- c('sample_id', 'cell_line', 'day_numeric')
  missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
  if (length(missing_metadata_columns) > 0L) {
    stop(
      'Missing validation metadata column(s): ',
      paste(missing_metadata_columns, collapse = ', '),
      '.',
      call. = FALSE
    )
  }
  sample_ids <- as.character(metadata$sample_id)
  if (
    anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids) ||
      !identical(colnames(counts), sample_ids) ||
      !identical(colnames(vst), sample_ids)
  ) {
    stop('Validation matrices and metadata must have identical, unique sample IDs.', call. = FALSE)
  }
  if (!identical(rownames(counts), rownames(vst))) {
    stop('Validation count and VST matrices must have identical gene rows.', call. = FALSE)
  }

  ordered_days <- sort(unique(metadata$day_numeric))
  heldout_cell_lines <- sort(unique(as.character(metadata$cell_line)))

  split_temporal_genes <- function(temporal_genes, training_metadata) {
    gene_day_means <- sapply(ordered_days, function(day_value) {
      day_samples <- training_metadata$sample_id[
        training_metadata$day_numeric == day_value
      ]
      rowMeans(vst[temporal_genes, day_samples, drop = FALSE], na.rm = TRUE)
    })
    gene_day_correlation <- apply(gene_day_means, 1L, function(gene_values) {
      stats::cor(
        ordered_days,
        gene_values,
        method = 'spearman',
        use = 'pairwise.complete.obs'
      )
    })
    gene_day_correlation[!is.finite(gene_day_correlation)] <- 0

    list(
      `All temporal` = temporal_genes,
      Maturation = names(gene_day_correlation)[gene_day_correlation > 0],
      Progenitor = names(gene_day_correlation)[gene_day_correlation < 0]
    )
  }

  score_gene_set <- function(
      gene_ids,
      gene_set_name,
      heldout_cell_line,
      training_metadata) {
    timing_fit <- score_differentiation_timing(
      expression_matrix = vst,
      metadata = metadata,
      temporal_genes = gene_ids,
      sample_id_col = 'sample_id',
      time_col = 'day_numeric',
      reference_samples = training_metadata$sample_id
    )
    heldout_scores <- timing_fit$scores[!timing_fit$scores$is_reference, , drop = FALSE]
    retained_pca_variance <- timing_fit$pca_variance[
      timing_fit$pca_variance$retained,
      ,
      drop = FALSE
    ]

    score_table <- data.frame(
      heldout_cell_line = heldout_cell_line,
      sample_id = heldout_scores$sample_id,
      actual_day = heldout_scores$observed_time,
      gene_set = gene_set_name,
      predicted_day = heldout_scores$predicted_time,
      residual = heldout_scores$predicted_time - heldout_scores$observed_time,
      stringsAsFactors = FALSE
    )

    list(
      gene_set = gene_set_name,
      gene_count = length(timing_fit$temporal_genes),
      n_pcs = timing_fit$n_pcs,
      retained_variance_percent = utils::tail(
        retained_pca_variance$cumulative_variance_percent,
        1L
      ),
      pca_variance_percent = retained_pca_variance$variance_percent,
      scores = score_table
    )
  }

  analyze_heldout_cell_line <- function(heldout_cell_line) {
    message('Analyzing held-out cell line: ', heldout_cell_line)
    training_metadata <- metadata[
      as.character(metadata$cell_line) != heldout_cell_line,
      ,
      drop = FALSE
    ]
    training_metadata$cell_line <- droplevels(factor(training_metadata$cell_line))
    training_sample_ids <- training_metadata$sample_id

    temporal_selection <- select_temporal_genes(
      raw_counts = counts[, training_sample_ids, drop = FALSE],
      vst_expression = vst[, training_sample_ids, drop = FALSE],
      metadata = training_metadata,
      sample_id_col = 'sample_id',
      time_col = 'day_numeric',
      adjustment_covariates = 'cell_line',
      expression_cpm_cutoff = expression_cpm_cutoff,
      lrt_padj_cutoff = lrt_padj_cutoff,
      vst_dynamic_range_cutoff = vst_dynamic_range_cutoff
    )
    gene_sets <- split_temporal_genes(
      temporal_genes = temporal_selection$temporal_genes,
      training_metadata = training_metadata
    )
    gene_set_results <- Map(
      score_gene_set,
      gene_sets,
      names(gene_sets),
      MoreArgs = list(
        heldout_cell_line = heldout_cell_line,
        training_metadata = training_metadata
      )
    )
    score_table <- do.call(rbind, lapply(gene_set_results, function(result) {
      result$scores
    }))
    pca_summary <- do.call(rbind, lapply(gene_set_results, function(result) {
      first_three_variance <- rep(NA_real_, 3L)
      available_components <- seq_len(min(3L, length(result$pca_variance_percent)))
      first_three_variance[available_components] <- result$pca_variance_percent[
        available_components
      ]
      data.frame(
        gene_set = result$gene_set,
        gene_count = result$gene_count,
        n_pcs = result$n_pcs,
        retained_variance_percent = round(result$retained_variance_percent, 2),
        pc1_percent = round(first_three_variance[[1L]], 2),
        pc2_percent = round(first_three_variance[[2L]], 2),
        pc3_percent = round(first_three_variance[[3L]], 2),
        stringsAsFactors = FALSE
      )
    }))

    list(
      heldout_cell_line = heldout_cell_line,
      temporal_genes = temporal_selection$temporal_genes,
      temporal_selection_summary = temporal_selection$summary,
      gene_sets = gene_sets,
      pca_summary = pca_summary,
      scores = score_table
    )
  }

  validation_results <- stats::setNames(
    lapply(heldout_cell_lines, analyze_heldout_cell_line),
    heldout_cell_lines
  )
  score_table <- do.call(rbind, lapply(validation_results, function(result) {
    result$scores
  }))
  rownames(score_table) <- NULL
  summary_table <- do.call(rbind, lapply(validation_results, function(result) {
    summary <- result$pca_summary
    summary$heldout_cell_line <- result$heldout_cell_line
    summary$n_temporal_genes <- length(result$temporal_genes)
    summary[
      ,
      c(
        'heldout_cell_line',
        'gene_set',
        'gene_count',
        'n_temporal_genes',
        'n_pcs',
        'retained_variance_percent',
        'pc1_percent',
        'pc2_percent',
        'pc3_percent'
      ),
      drop = FALSE
    ]
  }))
  rownames(summary_table) <- NULL

  list(
    settings = .leave_one_cell_line_out_validation_settings(
      expression_cpm_cutoff,
      lrt_padj_cutoff,
      vst_dynamic_range_cutoff
    ),
    scores = score_table,
    summary = summary_table,
    results = validation_results
  )
}
