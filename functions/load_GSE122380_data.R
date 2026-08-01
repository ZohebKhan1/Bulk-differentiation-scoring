# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/load_GSE122380_data.R
#
# input data:
# - data/GSE122380_metadata.rds
# - data/GSE122380_counts.rds
# - data/GSE122380_vst.rds
# ----

# Load and align the three tracked tutorial inputs at their external-ingestion boundary
load_GSE122380_data <- function(
    metadata_path = 'data/GSE122380_metadata.rds',
    counts_path = 'data/GSE122380_counts.rds',
    vst_path = 'data/GSE122380_vst.rds') {
  metadata <- readRDS(metadata_path)
  counts <- readRDS(counts_path)
  vst <- readRDS(vst_path)

  if (!is.data.frame(metadata)) {
    stop('GSE122380 metadata must be a data frame.', call. = FALSE)
  }
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop('GSE122380 counts must be a numeric matrix.', call. = FALSE)
  }
  if (!is.matrix(vst) || !is.numeric(vst)) {
    stop('GSE122380 VST expression must be a numeric matrix.', call. = FALSE)
  }

  required_metadata_columns <- c('sample_id', 'day_numeric', 'cell_line')
  missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
  if (length(missing_metadata_columns) > 0L) {
    stop(
      'GSE122380 metadata is missing columns: ',
      paste(missing_metadata_columns, collapse = ', '),
      '.',
      call. = FALSE
    )
  }

  sample_ids <- as.character(metadata$sample_id)
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop('GSE122380 metadata sample_id values must be nonblank and unique.', call. = FALSE)
  }
  if (
    !is.numeric(metadata$day_numeric) || anyNA(metadata$day_numeric) ||
      any(!is.finite(metadata$day_numeric))
  ) {
    stop('GSE122380 metadata day_numeric values must be finite numeric values.', call. = FALSE)
  }
  cell_lines <- as.character(metadata$cell_line)
  if (anyNA(cell_lines) || any(!nzchar(cell_lines))) {
    stop('GSE122380 metadata cell_line values must be nonblank.', call. = FALSE)
  }

  counts_sample_ids <- colnames(counts)
  vst_sample_ids <- colnames(vst)
  if (
    is.null(counts_sample_ids) || is.null(vst_sample_ids) ||
      anyNA(counts_sample_ids) || anyNA(vst_sample_ids) ||
      any(!nzchar(counts_sample_ids)) || any(!nzchar(vst_sample_ids)) ||
      anyDuplicated(counts_sample_ids) || anyDuplicated(vst_sample_ids)
  ) {
    stop('GSE122380 matrix column names must be nonblank, unique sample IDs.', call. = FALSE)
  }
  if (!setequal(sample_ids, counts_sample_ids) || !setequal(sample_ids, vst_sample_ids)) {
    stop('GSE122380 metadata, counts, and VST objects must contain the same samples.', call. = FALSE)
  }

  counts_gene_ids <- rownames(counts)
  vst_gene_ids <- rownames(vst)
  if (
    is.null(counts_gene_ids) || is.null(vst_gene_ids) ||
      anyNA(counts_gene_ids) || anyNA(vst_gene_ids) ||
      any(!nzchar(counts_gene_ids)) || any(!nzchar(vst_gene_ids)) ||
      anyDuplicated(counts_gene_ids) || anyDuplicated(vst_gene_ids)
  ) {
    stop('GSE122380 matrix row names must be nonblank, unique gene IDs.', call. = FALSE)
  }
  if (!identical(counts_gene_ids, vst_gene_ids)) {
    stop('GSE122380 counts and VST matrices must have identical ordered gene IDs.', call. = FALSE)
  }
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
    stop('GSE122380 counts must be finite, nonnegative integers.', call. = FALSE)
  }
  if (anyNA(vst) || any(!is.finite(vst))) {
    stop('GSE122380 VST values must be finite and complete.', call. = FALSE)
  }

  counts <- counts[, sample_ids, drop = FALSE]
  vst <- vst[, sample_ids, drop = FALSE]
  metadata$sample_id <- sample_ids
  metadata$day_factor <- factor(
    metadata$day_numeric,
    levels = sort(unique(metadata$day_numeric))
  )
  metadata$cell_line <- droplevels(factor(metadata$cell_line))

  input_paths <- c(
    metadata = metadata_path,
    counts = counts_path,
    vst = vst_path
  )

  input_md5 <- unname(tools::md5sum(input_paths))
  names(input_md5) <- names(input_paths)

  list(
    metadata = metadata,
    counts = counts,
    vst = vst,
    input_paths = input_paths,
    input_md5 = input_md5
  )
}
