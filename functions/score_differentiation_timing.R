# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/score_differentiation_timing.R
# ----

#' Score samples against an ordered reference differentiation trajectory
#'
#' `expression_matrix` must contain genes in rows and uniquely named samples in
#' columns. `metadata` must contain one unique row per expression sample.
#' `temporal_genes` are matched to unique gene row names; genes absent from the
#' matrix are omitted and at least two genes must remain.
#'
#' PCA is trained on the reference samples after gene centering without scaling.
#' Reference timepoint centroids are joined in ascending order. Each sample is
#' projected to its nearest finite polyline segment, so predicted times and
#' normalized scores are bounded by the first and last reference timepoints.
#'
#' The returned list contains `scores`, `pca_coordinates`, `pca_fit`,
#' `centroid_polyline`, `temporal_genes`, and `reference_time_range`. The score
#' table always uses `sample_id`, `observed_time`, `predicted_time`, and
#' `differentiation_score` regardless of the input metadata column names.
#'
#' The function has no persistence, console, random-state, or working-directory
#' side effects.
score_differentiation_timing <- function(
    expression_matrix,
    metadata,
    temporal_genes,
    sample_id_col = 'sample_id',
    time_col = 'day_numeric',
    reference_col = NULL,
    reference_values = NULL,
    n_pcs = 3L) {
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop('expression_matrix must be a numeric matrix or data frame.', call. = FALSE)
  }
  if (is.data.frame(expression_matrix)) {
    numeric_columns <- vapply(expression_matrix, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      stop('Every expression_matrix column must be numeric.', call. = FALSE)
    }
  } else if (!is.numeric(expression_matrix)) {
    stop('expression_matrix must be numeric.', call. = FALSE)
  }
  expression_matrix <- as.matrix(expression_matrix)

  gene_ids <- rownames(expression_matrix)
  sample_ids <- colnames(expression_matrix)
  if (
    is.null(gene_ids) || anyNA(gene_ids) || any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids)
  ) {
    stop('expression_matrix row names must be nonblank, unique gene IDs.', call. = FALSE)
  }
  if (
    is.null(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)
  ) {
    stop('expression_matrix column names must be nonblank, unique sample IDs.', call. = FALSE)
  }
  if (!is.data.frame(metadata)) {
    stop('metadata must be a data frame.', call. = FALSE)
  }

  required_columns <- c(sample_id_col, time_col)
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      'metadata is missing columns: ',
      paste(missing_columns, collapse = ', '),
      '.',
      call. = FALSE
    )
  }
  metadata_sample_ids <- as.character(metadata[[sample_id_col]])
  if (
    anyNA(metadata_sample_ids) || any(!nzchar(metadata_sample_ids)) ||
      anyDuplicated(metadata_sample_ids)
  ) {
    stop('metadata sample IDs must be nonblank and unique.', call. = FALSE)
  }
  metadata_rows <- match(sample_ids, metadata_sample_ids)
  if (anyNA(metadata_rows)) {
    missing_sample_ids <- sample_ids[is.na(metadata_rows)]
    stop(
      'metadata is missing expression samples: ',
      paste(utils::head(missing_sample_ids, 5L), collapse = ', '),
      '.',
      call. = FALSE
    )
  }
  metadata <- metadata[metadata_rows, , drop = FALSE]
  metadata[[sample_id_col]] <- metadata_sample_ids[metadata_rows]

  if (!is.numeric(metadata[[time_col]])) {
    stop('time_col must be numeric.', call. = FALSE)
  }
  if (!is.character(temporal_genes)) {
    stop('temporal_genes must be a character vector of gene IDs.', call. = FALSE)
  }
  if (anyNA(temporal_genes) || any(!nzchar(temporal_genes))) {
    stop('temporal_genes must contain nonblank gene IDs.', call. = FALSE)
  }
  retained_temporal_genes <- intersect(unique(temporal_genes), gene_ids)
  if (length(retained_temporal_genes) < 2L) {
    stop('At least two temporal genes must be present in expression_matrix.', call. = FALSE)
  }
  retained_expression <- expression_matrix[retained_temporal_genes, , drop = FALSE]
  if (anyNA(retained_expression) || any(!is.finite(retained_expression))) {
    stop('Expression values for retained temporal genes must be finite and complete.', call. = FALSE)
  }

  if (is.null(reference_col)) {
    reference_idx <- rep(TRUE, nrow(metadata))
  } else {
    if (!reference_col %in% names(metadata)) {
      stop('reference_col is not present in metadata.', call. = FALSE)
    }
    if (is.null(reference_values) || length(reference_values) == 0L) {
      stop('reference_values must be supplied when reference_col is used.', call. = FALSE)
    }
    reference_idx <- metadata[[reference_col]] %in% reference_values
  }
  if (sum(reference_idx) < 3L) {
    stop('At least three reference samples are required.', call. = FALSE)
  }
  reference_times <- metadata[[time_col]][reference_idx]
  if (anyNA(reference_times) || any(!is.finite(reference_times))) {
    stop('Reference timepoints must be finite and complete.', call. = FALSE)
  }
  reference_timepoints <- sort(unique(reference_times))
  if (length(reference_timepoints) < 2L) {
    stop('Reference samples must span at least two timepoints.', call. = FALSE)
  }

  if (
    length(n_pcs) != 1L || !is.numeric(n_pcs) || is.na(n_pcs) ||
      !is.finite(n_pcs) || n_pcs != as.integer(n_pcs)
  ) {
    stop('n_pcs must be one finite integer.', call. = FALSE)
  }
  n_pcs <- as.integer(n_pcs)
  max_pcs <- min(sum(reference_idx) - 1L, length(retained_temporal_genes))
  if (n_pcs < 1L || n_pcs > max_pcs) {
    stop('n_pcs must be between 1 and ', max_pcs, ' for these inputs.', call. = FALSE)
  }

  reference_matrix <- t(retained_expression[, reference_idx, drop = FALSE])
  pca_fit <- stats::prcomp(reference_matrix, center = TRUE, scale. = FALSE)
  rotation <- pca_fit$rotation[, seq_len(n_pcs), drop = FALSE]
  centered_matrix <- scale(
    t(retained_expression),
    center = pca_fit$center,
    scale = FALSE
  )
  pca_coordinates <- as.data.frame(centered_matrix %*% rotation)
  pc_columns <- paste0('PC', seq_len(n_pcs))
  colnames(pca_coordinates) <- pc_columns
  pca_coordinates[[sample_id_col]] <- rownames(pca_coordinates)
  pca_coordinates[[time_col]] <- metadata[[time_col]]
  pca_coordinates$is_reference <- reference_idx

  centroid_polyline <- stats::aggregate(
    pca_coordinates[reference_idx, pc_columns, drop = FALSE],
    by = list(timepoint = reference_times),
    FUN = mean
  )
  centroid_polyline <- centroid_polyline[
    order(centroid_polyline$timepoint),
    ,
    drop = FALSE
  ]
  projection <- .score_timing_project_to_polyline(
    point_matrix = as.matrix(pca_coordinates[, pc_columns, drop = FALSE]),
    centroid_polyline = centroid_polyline,
    pc_columns = pc_columns
  )

  first_timepoint <- centroid_polyline$timepoint[[1]]
  last_timepoint <- centroid_polyline$timepoint[[nrow(centroid_polyline)]]
  time_span <- last_timepoint - first_timepoint

  scores <- data.frame(
    sample_id = metadata[[sample_id_col]],
    observed_time = metadata[[time_col]],
    is_reference = reference_idx,
    predicted_time = projection$predicted_time,
    differentiation_score = (projection$predicted_time - first_timepoint) / time_span,
    nearest_segment_start = projection$segment_start,
    nearest_segment_end = projection$segment_end,
    segment_fraction = projection$segment_fraction,
    squared_distance = projection$squared_distance,
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    pca_coordinates = pca_coordinates,
    pca_fit = pca_fit,
    centroid_polyline = centroid_polyline,
    temporal_genes = retained_temporal_genes,
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
    stop('The reference centroid polyline has no nonzero-length segments.', call. = FALSE)
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
