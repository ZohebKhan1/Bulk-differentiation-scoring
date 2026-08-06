#' Score samples against an ordered reference differentiation trajectory
#'
#' `expression_matrix` must contain genes in rows and uniquely named samples in
#' columns. `metadata` must contain one unique row per expression sample.
#' `temporal_genes` are matched to gene row names; absent genes are omitted.
#'
#' PCA is trained on the reference samples after gene centering without scaling.
#' Reference timepoint centroids are joined in ascending order. Each sample is
#' projected to its nearest finite polyline segment, so predicted times and
#' normalized scores are bounded by the first and last reference timepoints.
#'
#' @param expression_matrix Numeric genes-by-samples expression matrix.
#' @param metadata Data frame containing sample IDs and numeric timepoints.
#' @param temporal_genes Character vector of temporal-gene IDs.
#' @param sample_id_col Metadata column containing unique sample IDs.
#' @param time_col Metadata column containing finite numeric timepoints.
#' @param reference_col Optional metadata column identifying reference samples.
#' @param reference_values Values in `reference_col` that define the reference.
#' @param n_pcs Number of centered, unscaled PCA dimensions used for projection.
#'
#' @return A list containing `scores`, `pca_coordinates`, `pca_fit`,
#'   `centroid_polyline`, retained `temporal_genes`, and `reference_time_range`.
#'   The score table uses stable output names regardless of input column names.
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
  expression_matrix <- as.matrix(expression_matrix)

  gene_ids <- rownames(expression_matrix)
  sample_ids <- colnames(expression_matrix)
  metadata_sample_ids <- as.character(metadata[[sample_id_col]])
  metadata_rows <- match(sample_ids, metadata_sample_ids)
  metadata <- metadata[metadata_rows, , drop = FALSE]
  metadata[[sample_id_col]] <- metadata_sample_ids[metadata_rows]

  sample_times <- metadata[[time_col]]
  retained_temporal_genes <- intersect(unique(temporal_genes), gene_ids)
  retained_expression <- expression_matrix[retained_temporal_genes, , drop = FALSE]

  if (is.null(reference_col)) {
    reference_idx <- rep(TRUE, nrow(metadata))
  } else {
    reference_idx <- metadata[[reference_col]] %in% reference_values
  }
  reference_times <- sample_times[reference_idx]
  n_pcs <- as.integer(n_pcs)

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
