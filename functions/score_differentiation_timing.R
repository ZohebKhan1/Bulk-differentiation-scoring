#' Score samples against an ordered reference differentiation polyline
#'
#' `expression_matrix` must contain normalized expression with genes in rows and
#' uniquely named samples in columns. `metadata` must contain exactly one row per
#' expression sample. Temporal genes absent from the expression matrix are
#' omitted and reported in the returned object.
#'
#' The PCA is trained on the mean expression profile at each reference
#' timepoint. Genes are centered without variance scaling. The scorer retains
#' the smallest number of PCs explaining at least 99% of the between-timepoint
#' variance, joins the ordered reference-timepoint centroids with finite line
#' segments, and projects every sample to its nearest point on that polyline.
#'
#' @param expression_matrix Numeric genes-by-samples normalized expression
#'   matrix, such as variance-stabilized expression.
#' @param metadata Data frame containing sample IDs and numeric timepoints.
#'   Non-reference samples may have missing timepoints.
#' @param temporal_genes Character vector of temporal-gene IDs selected without
#'   using the samples being evaluated.
#' @param sample_id_col Metadata column containing unique sample IDs.
#' @param time_col Metadata column containing numeric reference timepoints.
#' @param reference_samples Optional character vector of sample IDs used to
#'   train the PCA and reference polyline. All samples are references when this
#'   is `NULL`.
#'
#' @return A list containing `scores`, individual-sample `pca_coordinates`, the
#'   full `pca_fit`, `pca_variance`, automatically selected `n_pcs`, the
#'   `centroid_polyline`, retained and omitted temporal genes, reference sample
#'   counts, and the reference time range. Output column names are stable
#'   regardless of the input metadata column names.
#'
#' The function has no persistence, console, random-state, or working-directory
#' side effects.
score_differentiation_timing <- function(
    expression_matrix,
    metadata,
    temporal_genes,
    sample_id_col = 'sample_id',
    time_col = 'day_numeric',
    reference_samples = NULL) {
  pca_variance_target <- 0.99

  expression_matrix <- as.matrix(expression_matrix)
  if (!is.numeric(expression_matrix) || length(dim(expression_matrix)) != 2L) {
    stop('`expression_matrix` must be a numeric matrix.', call. = FALSE)
  }
  if (nrow(expression_matrix) < 1L || ncol(expression_matrix) < 1L) {
    stop('`expression_matrix` must contain at least one gene and sample.', call. = FALSE)
  }

  gene_ids <- rownames(expression_matrix)
  sample_ids <- colnames(expression_matrix)
  if (is.null(gene_ids) || anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop('Expression-matrix gene row names must be present, non-empty, and unique.', call. = FALSE)
  }
  if (
    is.null(sample_ids) ||
      anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)
  ) {
    stop('Expression-matrix sample column names must be present, non-empty, and unique.', call. = FALSE)
  }

  if (!is.data.frame(metadata)) {
    stop('`metadata` must be a data frame.', call. = FALSE)
  }
  if (!is.character(sample_id_col) || length(sample_id_col) != 1L || !nzchar(sample_id_col)) {
    stop('`sample_id_col` must name one metadata column.', call. = FALSE)
  }
  if (!is.character(time_col) || length(time_col) != 1L || !nzchar(time_col)) {
    stop('`time_col` must name one metadata column.', call. = FALSE)
  }
  missing_metadata_columns <- setdiff(c(sample_id_col, time_col), names(metadata))
  if (length(missing_metadata_columns) > 0L) {
    stop(
      'Missing metadata column(s): ',
      paste(missing_metadata_columns, collapse = ', '),
      '.',
      call. = FALSE
    )
  }

  metadata_sample_ids <- as.character(metadata[[sample_id_col]])
  if (
    anyNA(metadata_sample_ids) ||
      any(!nzchar(metadata_sample_ids)) ||
      anyDuplicated(metadata_sample_ids)
  ) {
    stop('Metadata sample IDs must be present, non-empty, and unique.', call. = FALSE)
  }
  missing_metadata_samples <- setdiff(sample_ids, metadata_sample_ids)
  extra_metadata_samples <- setdiff(metadata_sample_ids, sample_ids)
  if (length(missing_metadata_samples) > 0L || length(extra_metadata_samples) > 0L) {
    stop(
      '`metadata` sample IDs must match expression-matrix column names exactly.',
      call. = FALSE
    )
  }
  metadata <- metadata[match(sample_ids, metadata_sample_ids), , drop = FALSE]

  sample_times <- metadata[[time_col]]
  if (!is.numeric(sample_times)) {
    stop('The metadata time column must be numeric.', call. = FALSE)
  }
  if (any(!is.na(sample_times) & !is.finite(sample_times))) {
    stop('The metadata time column cannot contain infinite values.', call. = FALSE)
  }

  if (!is.character(temporal_genes) || length(temporal_genes) < 1L) {
    stop('`temporal_genes` must be a non-empty character vector.', call. = FALSE)
  }
  temporal_genes <- unique(temporal_genes[!is.na(temporal_genes) & nzchar(temporal_genes)])
  retained_temporal_genes <- intersect(temporal_genes, gene_ids)
  missing_temporal_genes <- setdiff(temporal_genes, gene_ids)
  if (length(retained_temporal_genes) < 1L) {
    stop('No temporal genes are present in `expression_matrix`.', call. = FALSE)
  }

  if (is.null(reference_samples)) {
    reference_samples <- sample_ids
  } else {
    if (!is.atomic(reference_samples) || length(reference_samples) < 1L) {
      stop('`reference_samples` must be a non-empty vector of sample IDs.', call. = FALSE)
    }
    reference_samples <- unique(as.character(reference_samples))
    if (anyNA(reference_samples) || any(!nzchar(reference_samples))) {
      stop('`reference_samples` cannot contain missing or empty sample IDs.', call. = FALSE)
    }
    unknown_reference_samples <- setdiff(reference_samples, sample_ids)
    if (length(unknown_reference_samples) > 0L) {
      stop(
        'Unknown reference sample ID(s): ',
        paste(utils::head(unknown_reference_samples, 5L), collapse = ', '),
        if (length(unknown_reference_samples) > 5L) ', ...' else '',
        '.',
        call. = FALSE
      )
    }
  }
  reference_idx <- sample_ids %in% reference_samples
  reference_times <- sample_times[reference_idx]
  if (anyNA(reference_times) || any(!is.finite(reference_times))) {
    stop('Every reference sample must have a finite numeric timepoint.', call. = FALSE)
  }
  reference_timepoints <- sort(unique(reference_times))
  if (length(reference_timepoints) < 2L) {
    stop('At least two distinct reference timepoints are required.', call. = FALSE)
  }

  retained_expression <- expression_matrix[retained_temporal_genes, , drop = FALSE]
  if (any(!is.finite(retained_expression))) {
    stop('Retained temporal-gene expression values must all be finite.', call. = FALSE)
  }

  reference_day_means <- vapply(reference_timepoints, function(timepoint) {
    rowMeans(
      retained_expression[, reference_idx & sample_times == timepoint, drop = FALSE]
    )
  }, numeric(length(retained_temporal_genes)))
  rownames(reference_day_means) <- retained_temporal_genes
  reference_day_means <- t(reference_day_means)
  rownames(reference_day_means) <- as.character(reference_timepoints)

  between_timepoint_variance <- apply(reference_day_means, 2L, stats::var)
  informative_genes <- names(between_timepoint_variance)[
    is.finite(between_timepoint_variance) & between_timepoint_variance > 0
  ]
  invariant_temporal_genes <- setdiff(retained_temporal_genes, informative_genes)
  if (length(informative_genes) < 1L) {
    stop('Temporal genes have no variation between reference timepoint means.', call. = FALSE)
  }
  retained_temporal_genes <- informative_genes
  retained_expression <- retained_expression[retained_temporal_genes, , drop = FALSE]
  reference_day_means <- reference_day_means[, retained_temporal_genes, drop = FALSE]

  pca_fit <- stats::prcomp(reference_day_means, center = TRUE, scale. = FALSE)
  component_variance <- pca_fit$sdev^2
  total_variance <- sum(component_variance)
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop('Reference timepoint means do not define a finite PCA space.', call. = FALSE)
  }
  variance_proportion <- component_variance / total_variance
  cumulative_variance <- cumsum(variance_proportion)
  n_pcs <- which(cumulative_variance >= pca_variance_target)[[1L]]
  pc_columns <- paste0('PC', seq_len(n_pcs))
  rotation <- pca_fit$rotation[, pc_columns, drop = FALSE]

  centered_expression <- sweep(
    t(retained_expression),
    MARGIN = 2L,
    STATS = pca_fit$center,
    FUN = '-'
  )
  projected_expression <- centered_expression %*% rotation
  colnames(projected_expression) <- pc_columns

  centered_day_means <- sweep(
    reference_day_means,
    MARGIN = 2L,
    STATS = pca_fit$center,
    FUN = '-'
  )
  centroid_coordinates <- centered_day_means %*% rotation
  colnames(centroid_coordinates) <- pc_columns
  centroid_polyline <- data.frame(
    timepoint = reference_timepoints,
    centroid_coordinates,
    check.names = FALSE,
    row.names = NULL
  )

  projection <- .score_timing_project_to_polyline(
    point_matrix = projected_expression,
    centroid_polyline = centroid_polyline,
    pc_columns = pc_columns
  )

  first_timepoint <- reference_timepoints[[1L]]
  last_timepoint <- reference_timepoints[[length(reference_timepoints)]]
  time_span <- last_timepoint - first_timepoint
  scores <- data.frame(
    sample_id = sample_ids,
    observed_time = sample_times,
    is_reference = reference_idx,
    predicted_time = projection$predicted_time,
    differentiation_score = (projection$predicted_time - first_timepoint) / time_span,
    nearest_segment_start = projection$segment_start,
    nearest_segment_end = projection$segment_end,
    segment_fraction = projection$segment_fraction,
    squared_distance = projection$squared_distance,
    stringsAsFactors = FALSE
  )
  pca_coordinates <- data.frame(
    sample_id = sample_ids,
    observed_time = sample_times,
    is_reference = reference_idx,
    projected_expression,
    check.names = FALSE,
    row.names = NULL
  )
  pca_variance <- data.frame(
    component = paste0('PC', seq_along(component_variance)),
    variance_percent = variance_proportion * 100,
    cumulative_variance_percent = cumulative_variance * 100,
    retained = seq_along(component_variance) <= n_pcs,
    stringsAsFactors = FALSE
  )
  reference_timepoint_counts <- data.frame(
    timepoint = reference_timepoints,
    n_reference_samples = as.integer(
      table(factor(reference_times, levels = reference_timepoints))
    ),
    row.names = NULL
  )

  list(
    scores = scores,
    pca_coordinates = pca_coordinates,
    pca_fit = pca_fit,
    pca_variance = pca_variance,
    n_pcs = n_pcs,
    pca_variance_target = pca_variance_target,
    centroid_polyline = centroid_polyline,
    temporal_genes = retained_temporal_genes,
    missing_temporal_genes = missing_temporal_genes,
    invariant_temporal_genes = invariant_temporal_genes,
    reference_samples = sample_ids[reference_idx],
    reference_timepoint_counts = reference_timepoint_counts,
    reference_time_range = c(start = first_timepoint, end = last_timepoint)
  )
}

.score_timing_project_to_polyline <- function(point_matrix, centroid_polyline, pc_columns) {
  start_points <- as.matrix(
    centroid_polyline[-nrow(centroid_polyline), pc_columns, drop = FALSE]
  )
  end_points <- as.matrix(centroid_polyline[-1L, pc_columns, drop = FALSE])
  segment_vectors <- end_points - start_points
  segment_lengths_squared <- rowSums(segment_vectors^2)
  keep_segments <- is.finite(segment_lengths_squared) & segment_lengths_squared > 0
  if (!any(keep_segments)) {
    stop('The reference centroids do not define a non-zero polyline segment.', call. = FALSE)
  }

  start_points <- start_points[keep_segments, , drop = FALSE]
  segment_vectors <- segment_vectors[keep_segments, , drop = FALSE]
  segment_lengths_squared <- segment_lengths_squared[keep_segments]
  segment_start <- centroid_polyline$timepoint[-nrow(centroid_polyline)][keep_segments]
  segment_end <- centroid_polyline$timepoint[-1L][keep_segments]

  projected <- lapply(seq_len(nrow(point_matrix)), function(point_index) {
    point <- point_matrix[point_index, ]
    point_repeated <- matrix(
      point,
      nrow = nrow(start_points),
      ncol = ncol(start_points),
      byrow = TRUE
    )
    segment_offset <- point_repeated - start_points
    segment_fraction <- rowSums(segment_offset * segment_vectors) / segment_lengths_squared
    segment_fraction <- pmin(pmax(segment_fraction, 0), 1)
    projected_points <- start_points + segment_fraction * segment_vectors
    squared_distance <- rowSums((point_repeated - projected_points)^2)
    best_segment <- which.min(squared_distance)
    predicted_time <- segment_start[[best_segment]] +
      segment_fraction[[best_segment]] *
        (segment_end[[best_segment]] - segment_start[[best_segment]])

    data.frame(
      predicted_time = predicted_time,
      segment_start = segment_start[[best_segment]],
      segment_end = segment_end[[best_segment]],
      segment_fraction = segment_fraction[[best_segment]],
      squared_distance = squared_distance[[best_segment]]
    )
  })

  do.call(rbind, projected)
}
