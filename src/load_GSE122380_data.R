#' Load and align the tracked GSE122380 inputs
#'
#' @return A list containing aligned `metadata`, `counts`, and `vst` objects.
load_GSE122380_data <- function(
    metadata_path = 'data/GSE122380_metadata.rds',
    counts_path = 'data/GSE122380_counts.rds',
    vst_path = 'data/GSE122380_vst.rds') {
  metadata <- readRDS(metadata_path)
  counts <- readRDS(counts_path)
  vst <- readRDS(vst_path)

  sample_ids <- as.character(metadata$sample_id)
  counts <- counts[, sample_ids, drop = FALSE]
  vst <- vst[, sample_ids, drop = FALSE]
  metadata$sample_id <- sample_ids
  metadata$day_factor <- factor(
    metadata$day_numeric,
    levels = sort(unique(metadata$day_numeric))
  )
  metadata$cell_line <- droplevels(factor(metadata$cell_line))

  list(
    metadata = metadata,
    counts = counts,
    vst = vst
  )
}
