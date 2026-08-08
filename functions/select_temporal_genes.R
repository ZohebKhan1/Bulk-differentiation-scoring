#' Select temporal genes from raw counts and VST expression
#'
#' Genes must pass a maximum day-mean TMM CPM threshold, a DESeq2 likelihood-
#' ratio test (LRT) for categorical time, and a minimum day-mean VST range.
#' The LRT full model contains the requested adjustment covariates plus time;
#' the reduced model contains the adjustment covariates alone. An optional
#' metadata column/value pair restricts all three filters to a reference group.
#'
#' @param raw_counts Integer-like raw-count matrix with genes in rows and
#'   samples in columns. TMM CPM is calculated internally from these counts;
#'   normalized CPM values must not be supplied here.
#' @param vst_expression Numeric variance-stabilized expression matrix with the
#'   same genes and samples as `raw_counts`.
#' @param metadata Sample-level data frame containing one row per matrix sample.
#' @param sample_id_col Metadata column containing unique sample identifiers.
#' @param time_col Metadata column containing finite numeric timepoints. Time is
#'   treated as a categorical factor in the DESeq2 LRT.
#' @param adjustment_covariates Character vector of metadata columns included
#'   in both the full and reduced DESeq2 models. Use `NULL` for an unadjusted
#'   time LRT, which is the default.
#' @param reference_group_col Optional metadata column used to select the
#'   reference cohort. Supply this together with `reference_group_value`, or
#'   leave both as `NULL` to use all samples.
#' @param reference_group_value Single value in `reference_group_col` that
#'   identifies reference samples.
#' @param expression_cpm_cutoff Minimum maximum day-mean TMM CPM.
#' @param lrt_padj_cutoff Maximum adjusted p-value from the DESeq2 LRT.
#' @param vst_dynamic_range_cutoff Minimum range across day-mean VST values.
#'
#' @return A list containing ordered temporal genes, LRT results, day-mean TMM
#'   CPM and VST matrices, the reference samples, the generated DESeq2 design,
#'   and a one-row filtering summary.
select_temporal_genes <- function(
  raw_counts,
  vst_expression,
  metadata,
  sample_id_col = 'sample_id',
  time_col = 'day_numeric',
  adjustment_covariates = NULL,
  reference_group_col = NULL,
  reference_group_value = NULL,
  expression_cpm_cutoff = 10,
  lrt_padj_cutoff = 1e-7,
  vst_dynamic_range_cutoff = 0.6) {
  raw_counts <- .select_temporal_validate_matrix(raw_counts, 'raw_counts')
  vst_expression <- .select_temporal_validate_matrix(
    vst_expression,
    'vst_expression'
  )

  if (any(raw_counts < 0)) {
    stop('`raw_counts` cannot contain negative values.', call. = FALSE)
  }
  if (any(abs(raw_counts - round(raw_counts)) > 1e-8)) {
    stop(
      '`raw_counts` must contain integer-like, unnormalized counts; do not supply CPM values.',
      call. = FALSE
    )
  }

  .select_temporal_validate_threshold(
    expression_cpm_cutoff,
    'expression_cpm_cutoff',
    lower_bound = 0,
    lower_inclusive = TRUE
  )
  .select_temporal_validate_threshold(
    lrt_padj_cutoff,
    'lrt_padj_cutoff',
    lower_bound = 0,
    upper_bound = 1,
    lower_inclusive = FALSE,
    upper_inclusive = FALSE
  )
  .select_temporal_validate_threshold(
    vst_dynamic_range_cutoff,
    'vst_dynamic_range_cutoff',
    lower_bound = 0,
    lower_inclusive = TRUE
  )

  .select_temporal_validate_column_name(sample_id_col, 'sample_id_col')
  .select_temporal_validate_column_name(time_col, 'time_col')
  if (!is.data.frame(metadata) || nrow(metadata) < 1L) {
    stop('`metadata` must be a non-empty data frame.', call. = FALSE)
  }
  if (anyDuplicated(names(metadata))) {
    stop('`metadata` column names must be unique.', call. = FALSE)
  }

  if (is.null(adjustment_covariates)) {
    adjustment_covariates <- character()
  }
  if (
    !is.character(adjustment_covariates) ||
      anyNA(adjustment_covariates) ||
      any(!nzchar(adjustment_covariates))
  ) {
    stop(
      '`adjustment_covariates` must be `NULL` or a character vector of metadata column names.',
      call. = FALSE
    )
  }
  if (anyDuplicated(adjustment_covariates)) {
    stop('`adjustment_covariates` cannot contain duplicate column names.', call. = FALSE)
  }
  forbidden_covariates <- intersect(
    adjustment_covariates,
    c(sample_id_col, time_col)
  )
  if (length(forbidden_covariates) > 0L) {
    stop(
      '`adjustment_covariates` cannot include the sample-ID or time column: ',
      .select_temporal_format_values(forbidden_covariates),
      '. Time is added to the full model automatically.',
      call. = FALSE
    )
  }

  reference_arguments_supplied <- c(
    !is.null(reference_group_col),
    !is.null(reference_group_value)
  )
  if (sum(reference_arguments_supplied) == 1L) {
    stop(
      '`reference_group_col` and `reference_group_value` must be supplied together, or both left as `NULL`.',
      call. = FALSE
    )
  }
  if (!is.null(reference_group_col)) {
    .select_temporal_validate_column_name(
      reference_group_col,
      'reference_group_col'
    )
    if (
      length(reference_group_value) != 1L ||
        !is.atomic(reference_group_value) ||
        is.na(reference_group_value) ||
        is.numeric(reference_group_value) && !is.finite(reference_group_value) ||
        is.character(reference_group_value) && !nzchar(reference_group_value)
    ) {
      stop(
        '`reference_group_value` must be one non-missing atomic value.',
        call. = FALSE
      )
    }
    if (reference_group_col %in% adjustment_covariates) {
      stop(
        '`reference_group_col` cannot also be an adjustment covariate because filtering to one reference value makes it constant.',
        call. = FALSE
      )
    }
  }

  missing_core_columns <- setdiff(c(sample_id_col, time_col), names(metadata))
  if (length(missing_core_columns) > 0L) {
    stop(
      '`metadata` is missing the requested sample-ID or time column(s): ',
      .select_temporal_format_values(missing_core_columns),
      '. Available columns: ',
      .select_temporal_format_values(names(metadata), max_values = 12L),
      '.',
      call. = FALSE
    )
  }
  if (
    !is.null(reference_group_col) &&
      !reference_group_col %in% names(metadata)
  ) {
    stop(
      '`reference_group_col` specifies "',
      reference_group_col,
      '", but `metadata` does not contain that column. Available columns: ',
      .select_temporal_format_values(names(metadata), max_values = 12L),
      '.',
      call. = FALSE
    )
  }
  missing_adjustment_covariates <- setdiff(
    adjustment_covariates,
    names(metadata)
  )
  if (length(missing_adjustment_covariates) > 0L) {
    stop(
      '`adjustment_covariates` contains column(s) not found in `metadata`: ',
      .select_temporal_format_values(missing_adjustment_covariates),
      '. Available columns: ',
      .select_temporal_format_values(names(metadata), max_values = 12L),
      '. Set `adjustment_covariates = NULL` for an unadjusted time LRT.',
      call. = FALSE
    )
  }

  metadata_sample_ids <- metadata[[sample_id_col]]
  if (!is.atomic(metadata_sample_ids) && !is.factor(metadata_sample_ids)) {
    stop(
      '`metadata[["',
      sample_id_col,
      '"]]` must be an atomic vector or factor.',
      call. = FALSE
    )
  }
  sample_ids <- as.character(metadata_sample_ids)
  invalid_sample_ids <- is.na(sample_ids) | !nzchar(sample_ids)
  if (any(invalid_sample_ids)) {
    stop(
      '`metadata[["',
      sample_id_col,
      '"]]` contains missing or empty sample IDs at row(s): ',
      .select_temporal_format_values(which(invalid_sample_ids)),
      '.',
      call. = FALSE
    )
  }
  duplicate_sample_ids <- unique(sample_ids[duplicated(sample_ids)])
  if (length(duplicate_sample_ids) > 0L) {
    stop(
      '`metadata[["',
      sample_id_col,
      '"]]` contains duplicate sample ID(s): ',
      .select_temporal_format_values(duplicate_sample_ids),
      '.',
      call. = FALSE
    )
  }

  raw_gene_ids <- rownames(raw_counts)
  vst_gene_ids <- rownames(vst_expression)
  missing_vst_genes <- setdiff(raw_gene_ids, vst_gene_ids)
  extra_vst_genes <- setdiff(vst_gene_ids, raw_gene_ids)
  if (length(missing_vst_genes) > 0L || length(extra_vst_genes) > 0L) {
    stop(
      '`raw_counts` and `vst_expression` must contain exactly the same gene IDs. ',
      .select_temporal_set_difference_message(
        missing_vst_genes,
        extra_vst_genes,
        'missing from `vst_expression`',
        'present only in `vst_expression`'
      ),
      call. = FALSE
    )
  }

  raw_sample_ids <- colnames(raw_counts)
  vst_sample_ids <- colnames(vst_expression)
  missing_vst_samples <- setdiff(raw_sample_ids, vst_sample_ids)
  extra_vst_samples <- setdiff(vst_sample_ids, raw_sample_ids)
  if (length(missing_vst_samples) > 0L || length(extra_vst_samples) > 0L) {
    stop(
      '`raw_counts` and `vst_expression` must contain exactly the same sample IDs. ',
      .select_temporal_set_difference_message(
        missing_vst_samples,
        extra_vst_samples,
        'missing from `vst_expression`',
        'present only in `vst_expression`'
      ),
      call. = FALSE
    )
  }
  missing_metadata_samples <- setdiff(raw_sample_ids, sample_ids)
  extra_metadata_samples <- setdiff(sample_ids, raw_sample_ids)
  if (
    length(missing_metadata_samples) > 0L ||
      length(extra_metadata_samples) > 0L
  ) {
    stop(
      '`metadata` and the expression matrices must contain exactly the same sample IDs. ',
      .select_temporal_set_difference_message(
        missing_metadata_samples,
        extra_metadata_samples,
        'missing from `metadata`',
        'present only in `metadata`'
      ),
      call. = FALSE
    )
  }

  vst_expression <- vst_expression[
    raw_gene_ids,
    raw_sample_ids,
    drop = FALSE
  ]
  metadata <- metadata[match(raw_sample_ids, sample_ids), , drop = FALSE]
  metadata[[sample_id_col]] <- raw_sample_ids

  if (is.null(reference_group_col)) {
    reference_idx <- rep(TRUE, nrow(metadata))
  } else {
    reference_values <- metadata[[reference_group_col]]
    if (!is.atomic(reference_values) && !is.factor(reference_values)) {
      stop(
        '`metadata[["',
        reference_group_col,
        '"]]` must be an atomic vector or factor.',
        call. = FALSE
      )
    }
    invalid_reference_values <- is.na(reference_values)
    if (is.numeric(reference_values)) {
      invalid_reference_values <- invalid_reference_values |
        !is.finite(reference_values)
    }
    if (is.character(reference_values) || is.factor(reference_values)) {
      invalid_reference_values <- invalid_reference_values |
        !nzchar(as.character(reference_values))
    }
    if (any(invalid_reference_values)) {
      stop(
        '`metadata[["',
        reference_group_col,
        '"]]` contains missing or empty group labels for sample(s): ',
        .select_temporal_format_values(
          raw_sample_ids[invalid_reference_values]
        ),
        '.',
        call. = FALSE
      )
    }
    available_reference_values <- sort(unique(as.character(reference_values)))
    reference_idx <- as.character(reference_values) ==
      as.character(reference_group_value)
    if (!any(reference_idx)) {
      stop(
        '`reference_group_value` "',
        as.character(reference_group_value),
        '" was not found in `metadata[["',
        reference_group_col,
        '"]]`. Available values: ',
        .select_temporal_format_values(available_reference_values),
        '.',
        call. = FALSE
      )
    }
  }

  reference_metadata <- droplevels(metadata[reference_idx, , drop = FALSE])
  reference_sample_ids <- reference_metadata[[sample_id_col]]
  reference_raw_counts <- raw_counts[
    ,
    reference_sample_ids,
    drop = FALSE
  ]
  reference_vst_expression <- vst_expression[
    ,
    reference_sample_ids,
    drop = FALSE
  ]

  sample_times <- reference_metadata[[time_col]]
  if (!is.numeric(sample_times)) {
    stop(
      '`metadata[["',
      time_col,
      '"]]` must be numeric; time is converted to a categorical factor internally for the LRT.',
      call. = FALSE
    )
  }
  invalid_times <- is.na(sample_times) | !is.finite(sample_times)
  if (any(invalid_times)) {
    stop(
      '`metadata[["',
      time_col,
      '"]]` contains missing or non-finite timepoints for reference sample(s): ',
      .select_temporal_format_values(
        reference_sample_ids[invalid_times]
      ),
      '.',
      call. = FALSE
    )
  }
  days <- sort(unique(sample_times))
  if (length(days) < 2L) {
    stop(
      'Reference samples must span at least two distinct values in `metadata[["',
      time_col,
      '"]]`; found: ',
      .select_temporal_format_values(days),
      '.',
      call. = FALSE
    )
  }

  lrt_metadata <- reference_metadata
  rownames(lrt_metadata) <- reference_sample_ids
  for (covariate in adjustment_covariates) {
    covariate_values <- lrt_metadata[[covariate]]
    invalid_covariate_values <- is.na(covariate_values)
    if (is.character(covariate_values) || is.factor(covariate_values)) {
      invalid_covariate_values <- invalid_covariate_values |
        !nzchar(as.character(covariate_values))
    }
    if (is.numeric(covariate_values)) {
      invalid_covariate_values <- invalid_covariate_values |
        !is.finite(covariate_values)
    }
    if (any(invalid_covariate_values)) {
      stop(
        'Adjustment covariate `',
        covariate,
        '` contains missing, empty, or non-finite values for reference sample(s): ',
        .select_temporal_format_values(
          reference_sample_ids[invalid_covariate_values]
        ),
        '.',
        call. = FALSE
      )
    }
    valid_covariate_type <- is.factor(covariate_values) ||
      is.character(covariate_values) ||
      is.logical(covariate_values) ||
      is.numeric(covariate_values)
    if (!valid_covariate_type) {
      stop(
        'Adjustment covariate `',
        covariate,
        '` must be numeric, character, logical, or a factor.',
        call. = FALSE
      )
    }
    if (is.character(covariate_values) || is.logical(covariate_values)) {
      covariate_values <- factor(covariate_values)
    } else if (is.factor(covariate_values)) {
      covariate_values <- droplevels(covariate_values)
    }
    if (length(unique(covariate_values)) < 2L) {
      stop(
        'Adjustment covariate `',
        covariate,
        '` has only one value among the selected reference samples. Remove it from `adjustment_covariates` or choose a broader reference group.',
        call. = FALSE
      )
    }
    lrt_metadata[[covariate]] <- covariate_values
  }

  lrt_metadata[[time_col]] <- factor(sample_times, levels = days)
  full_design <- stats::reformulate(c(adjustment_covariates, time_col))
  reduced_design <- if (length(adjustment_covariates) > 0L) {
    stats::reformulate(adjustment_covariates)
  } else {
    stats::as.formula('~ 1')
  }

  .select_temporal_validate_design(full_design, lrt_metadata, 'full')
  .select_temporal_validate_design(reduced_design, lrt_metadata, 'reduced')

  library_sizes <- colSums(reference_raw_counts)
  if (any(!is.finite(library_sizes) | library_sizes <= 0)) {
    stop(
      '`raw_counts` has zero or non-finite library sizes for reference sample(s): ',
      .select_temporal_format_values(
        reference_sample_ids[!is.finite(library_sizes) | library_sizes <= 0]
      ),
      '.',
      call. = FALSE
    )
  }

  dge <- edgeR::DGEList(counts = reference_raw_counts)
  dge <- edgeR::normLibSizes(dge, method = 'TMM')
  tmm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  day_mean_tmm_cpm <- vapply(days, function(day_value) {
    day_samples <- reference_sample_ids[sample_times == day_value]
    rowMeans(tmm_cpm[, day_samples, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(reference_raw_counts)))
  rownames(day_mean_tmm_cpm) <- raw_gene_ids
  colnames(day_mean_tmm_cpm) <- paste0('D', days)
  max_day_mean_tmm_cpm <- apply(day_mean_tmm_cpm, 1L, max, na.rm = TRUE)
  expression_genes <- names(max_day_mean_tmm_cpm)[
    max_day_mean_tmm_cpm >= expression_cpm_cutoff
  ]
  if (length(expression_genes) < 1L) {
    stop(
      'No genes passed `expression_cpm_cutoff = ',
      format(expression_cpm_cutoff),
      '`. The largest observed maximum day-mean TMM CPM was ',
      format(max(max_day_mean_tmm_cpm), digits = 4L),
      '.',
      call. = FALSE
    )
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = reference_raw_counts[expression_genes, , drop = FALSE],
    colData = lrt_metadata,
    design = full_design
  )
  dds <- DESeq2::DESeq(
    dds,
    test = 'LRT',
    reduced = reduced_design,
    quiet = TRUE
  )

  lrt_results <- as.data.frame(
    DESeq2::results(dds, alpha = lrt_padj_cutoff)
  )
  lrt_results$gene_id <- rownames(lrt_results)
  lrt_results <- lrt_results[
    order(lrt_results$padj, lrt_results$pvalue, na.last = TRUE),
    ,
    drop = FALSE
  ]
  lrt_genes <- lrt_results$gene_id[
    !is.na(lrt_results$padj) & lrt_results$padj < lrt_padj_cutoff
  ]
  if (length(lrt_genes) < 1L) {
    stop(
      'No genes passed `lrt_padj_cutoff = ',
      format(lrt_padj_cutoff, scientific = TRUE),
      '` for the requested DESeq2 design. Review the reference group, design covariates, and threshold.',
      call. = FALSE
    )
  }

  day_mean_vst <- vapply(days, function(day_value) {
    day_samples <- reference_sample_ids[sample_times == day_value]
    rowMeans(
      reference_vst_expression[lrt_genes, day_samples, drop = FALSE],
      na.rm = TRUE
    )
  }, numeric(length(lrt_genes)))
  rownames(day_mean_vst) <- lrt_genes
  colnames(day_mean_vst) <- paste0('D', days)
  vst_dynamic_range <- apply(day_mean_vst, 1L, function(gene_values) {
    max(gene_values, na.rm = TRUE) - min(gene_values, na.rm = TRUE)
  })
  temporal_genes <- lrt_genes[
    vst_dynamic_range[lrt_genes] >= vst_dynamic_range_cutoff
  ]
  if (length(temporal_genes) < 1L) {
    stop(
      'No genes passed `vst_dynamic_range_cutoff = ',
      format(vst_dynamic_range_cutoff),
      '`. The largest observed day-mean VST range among LRT genes was ',
      format(max(vst_dynamic_range), digits = 4L),
      '.',
      call. = FALSE
    )
  }

  full_design_label <- .select_temporal_format_formula(full_design)
  reduced_design_label <- .select_temporal_format_formula(reduced_design)
  summary <- data.frame(
    input_genes = nrow(raw_counts),
    input_samples = ncol(raw_counts),
    reference_samples = length(reference_sample_ids),
    reference_timepoints = length(days),
    genes_passing_day_mean_tmm_cpm = length(expression_genes),
    genes_tested_by_lrt = nrow(lrt_results),
    genes_passing_lrt_padj = length(lrt_genes),
    genes_passing_all_filters = length(temporal_genes),
    expression_cpm_cutoff = expression_cpm_cutoff,
    lrt_padj_cutoff = lrt_padj_cutoff,
    vst_dynamic_range_cutoff = vst_dynamic_range_cutoff,
    full_design = full_design_label,
    reduced_design = reduced_design_label,
    reference_group_col = if (is.null(reference_group_col)) {
      NA_character_
    } else {
      reference_group_col
    },
    reference_group_value = if (is.null(reference_group_value)) {
      NA_character_
    } else {
      as.character(reference_group_value)
    },
    stringsAsFactors = FALSE
  )

  list(
    temporal_genes = temporal_genes,
    lrt_results = lrt_results,
    day_mean_tmm_cpm = day_mean_tmm_cpm,
    day_mean_vst = day_mean_vst,
    reference_samples = reference_sample_ids,
    design = list(
      sample_id_col = sample_id_col,
      time_col = time_col,
      adjustment_covariates = adjustment_covariates,
      full = full_design,
      reduced = reduced_design,
      reference_group_col = reference_group_col,
      reference_group_value = reference_group_value
    ),
    summary = summary
  )
}

.select_temporal_validate_matrix <- function(value, argument_name) {
  value <- as.matrix(value)
  if (!is.numeric(value) || length(dim(value)) != 2L) {
    stop('`', argument_name, '` must be a numeric matrix.', call. = FALSE)
  }
  if (nrow(value) < 1L || ncol(value) < 1L) {
    stop(
      '`',
      argument_name,
      '` must contain at least one gene and one sample.',
      call. = FALSE
    )
  }
  if (any(!is.finite(value))) {
    stop('`', argument_name, '` must contain only finite values.', call. = FALSE)
  }

  gene_ids <- rownames(value)
  sample_ids <- colnames(value)
  if (
    is.null(gene_ids) ||
      anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids)
  ) {
    stop(
      '`',
      argument_name,
      '` gene row names must be present, non-empty, and unique.',
      call. = FALSE
    )
  }
  if (
    is.null(sample_ids) ||
      anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)
  ) {
    stop(
      '`',
      argument_name,
      '` sample column names must be present, non-empty, and unique.',
      call. = FALSE
    )
  }
  value
}

.select_temporal_validate_column_name <- function(value, argument_name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(
      '`',
      argument_name,
      '` must name one metadata column.',
      call. = FALSE
    )
  }
}

.select_temporal_validate_threshold <- function(
  value,
  argument_name,
  lower_bound,
  upper_bound = Inf,
  lower_inclusive,
  upper_inclusive = TRUE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop('`', argument_name, '` must be one finite numeric value.', call. = FALSE)
  }
  lower_invalid <- if (lower_inclusive) {
    value < lower_bound
  } else {
    value <= lower_bound
  }
  upper_invalid <- if (upper_inclusive) {
    value > upper_bound
  } else {
    value >= upper_bound
  }
  if (lower_invalid || upper_invalid) {
    interval_start <- if (lower_inclusive) '[' else '('
    interval_end <- if (upper_inclusive) ']' else ')'
    stop(
      '`',
      argument_name,
      '` must be in ',
      interval_start,
      format(lower_bound),
      ', ',
      if (is.finite(upper_bound)) format(upper_bound) else 'Inf',
      interval_end,
      '.',
      call. = FALSE
    )
  }
}

.select_temporal_validate_design <- function(design, metadata, design_name) {
  model_matrix <- tryCatch(
    stats::model.matrix(design, data = metadata),
    error = function(error) {
      stop(
        'Could not construct the ',
        design_name,
        ' DESeq2 model matrix from ',
        .select_temporal_format_formula(design),
        ': ',
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (qr(model_matrix)$rank < ncol(model_matrix)) {
    stop(
      'The ',
      design_name,
      ' DESeq2 model matrix is not full rank for ',
      .select_temporal_format_formula(design),
      '. The selected adjustment covariates may be confounded with time or with each other.',
      call. = FALSE
    )
  }
  if (nrow(model_matrix) <= ncol(model_matrix)) {
    stop(
      'The ',
      design_name,
      ' DESeq2 model has ',
      ncol(model_matrix),
      ' coefficient(s) for ',
      nrow(model_matrix),
      ' reference sample(s), leaving no residual degrees of freedom. Use more reference samples or fewer adjustment covariates.',
      call. = FALSE
    )
  }
  invisible(model_matrix)
}

.select_temporal_format_formula <- function(formula) {
  paste(deparse(formula, width.cutoff = 500L), collapse = '')
}

.select_temporal_format_values <- function(values, max_values = 5L) {
  values <- as.character(values)
  displayed_values <- utils::head(values, max_values)
  paste0(
    paste0('"', displayed_values, '"', collapse = ', '),
    if (length(values) > max_values) ', ...' else ''
  )
}

.select_temporal_set_difference_message <- function(
  missing_values,
  extra_values,
  missing_label,
  extra_label) {
  messages <- character()
  if (length(missing_values) > 0L) {
    messages <- c(
      messages,
      paste0(
        missing_label,
        ': ',
        .select_temporal_format_values(missing_values)
      )
    )
  }
  if (length(extra_values) > 0L) {
    messages <- c(
      messages,
      paste0(
        extra_label,
        ': ',
        .select_temporal_format_values(extra_values)
      )
    )
  }
  paste(messages, collapse = '; ')
}
