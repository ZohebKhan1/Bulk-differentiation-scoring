source('functions/select_temporal_genes.R')

expect_error_message <- function(expression, expected_text) {
  observed_error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  if (is.null(observed_error)) {
    stop('Expected an error containing: ', expected_text, call. = FALSE)
  }
  observed_text <- conditionMessage(observed_error)
  if (!grepl(expected_text, observed_text, fixed = TRUE)) {
    stop(
      'Expected error containing "',
      expected_text,
      '", observed: ',
      observed_text,
      call. = FALSE
    )
  }
  invisible(observed_text)
}

stopifnot(
  identical(formals(select_temporal_genes)$expression_cpm_cutoff, 10),
  identical(formals(select_temporal_genes)$lrt_padj_cutoff, 1e-7),
  identical(formals(select_temporal_genes)$vst_dynamic_range_cutoff, 0.6)
)

set.seed(20260808)
metadata <- expand.grid(
  condition = c('control', 'treated'),
  cell_line = c('line_a', 'line_b', 'line_c'),
  day_numeric = 1:4,
  replicate = 1:2,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
metadata$sample_id <- sprintf('sample_%03d', seq_len(nrow(metadata)))
metadata <- metadata[
  ,
  c('sample_id', 'condition', 'cell_line', 'day_numeric', 'replicate')
]

n_genes <- 60L
gene_ids <- sprintf('gene_%03d', seq_len(n_genes))
gene_baseline <- 35 + seq_len(n_genes)
gene_time_slope <- c(seq(5, 24, length.out = 40L), rep(0, 20L))
cell_line_effect <- c(line_a = 0, line_b = 8, line_c = 16)
condition_effect <- c(control = 0, treated = 12)
expected_means <- vapply(seq_len(nrow(metadata)), function(sample_index) {
  gene_baseline +
    gene_time_slope * metadata$day_numeric[[sample_index]] +
    cell_line_effect[[metadata$cell_line[[sample_index]]]] +
    condition_effect[[metadata$condition[[sample_index]]]]
}, numeric(n_genes))
raw_counts <- matrix(
  stats::rnbinom(
    length(expected_means),
    mu = as.vector(expected_means),
    size = 30
  ),
  nrow = n_genes,
  dimnames = list(gene_ids, metadata$sample_id)
)
vst_expression <- log2(raw_counts + 1)

selection <- select_temporal_genes(
  raw_counts = raw_counts,
  vst_expression = vst_expression[rev(gene_ids), rev(metadata$sample_id)],
  metadata = metadata,
  adjustment_covariates = 'cell_line',
  reference_group_col = 'condition',
  reference_group_value = 'control',
  expression_cpm_cutoff = 0,
  lrt_padj_cutoff = 0.999,
  vst_dynamic_range_cutoff = 0
)

stopifnot(
  length(selection$temporal_genes) > 0L,
  length(selection$reference_samples) == sum(metadata$condition == 'control'),
  all(metadata$condition[
    match(selection$reference_samples, metadata$sample_id)
  ] == 'control'),
  identical(selection$design$adjustment_covariates, 'cell_line'),
  identical(selection$design$reference_group_col, 'condition'),
  identical(selection$design$reference_group_value, 'control'),
  selection$summary$reference_timepoints == 4L
)

expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    metadata[, c('sample_id', 'day_numeric')],
    reference_group_col = 'missing_condition',
    reference_group_value = 'control'
  ),
  '`reference_group_col` specifies "missing_condition", but `metadata` does not contain that column'
)

expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    metadata,
    reference_group_col = 'condition',
    reference_group_value = 'unknown'
  ),
  '`reference_group_value` "unknown" was not found'
)

expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    metadata,
    adjustment_covariates = 'missing_batch'
  ),
  '`adjustment_covariates` contains column(s) not found in `metadata`'
)

expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    metadata,
    reference_group_col = 'condition'
  ),
  '`reference_group_col` and `reference_group_value` must be supplied together'
)

constant_metadata <- metadata
constant_metadata$study <- 'only_study'
expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    constant_metadata,
    adjustment_covariates = 'study'
  ),
  'Adjustment covariate `study` has only one value'
)

missing_covariate_metadata <- metadata
missing_covariate_metadata$cell_line[[1L]] <- NA_character_
expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    missing_covariate_metadata
  ),
  'Adjustment covariate `cell_line` contains missing, empty, or non-finite values'
)

missing_time_metadata <- metadata
missing_time_metadata$day_numeric[[1L]] <- NA_real_
expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    missing_time_metadata
  ),
  '`metadata[["day_numeric"]]` contains missing or non-finite timepoints'
)

single_time_metadata <- metadata
single_time_metadata$day_numeric <- 1
expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    single_time_metadata
  ),
  'Reference samples must span at least two distinct values'
)

confounded_metadata <- metadata
confounded_metadata$batch <- paste0('day_', confounded_metadata$day_numeric)
expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    confounded_metadata,
    adjustment_covariates = 'batch'
  ),
  'DESeq2 model matrix is not full rank'
)

expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression[, -1L, drop = FALSE],
    metadata
  ),
  '`raw_counts` and `vst_expression` must contain exactly the same sample IDs'
)

duplicate_metadata <- metadata
duplicate_metadata$sample_id[[2L]] <- duplicate_metadata$sample_id[[1L]]
expect_error_message(
  select_temporal_genes(raw_counts, vst_expression, duplicate_metadata),
  '`metadata[["sample_id"]]` contains duplicate sample ID(s)'
)

noninteger_counts <- raw_counts
noninteger_counts[[1L]] <- noninteger_counts[[1L]] + 0.5
expect_error_message(
  select_temporal_genes(noninteger_counts, vst_expression, metadata),
  '`raw_counts` must contain integer-like, unnormalized counts; do not supply CPM values'
)

expect_error_message(
  select_temporal_genes(
    raw_counts,
    vst_expression,
    metadata,
    lrt_padj_cutoff = 0
  ),
  '`lrt_padj_cutoff` must be in (0, 1)'
)

message('select_temporal_genes tests passed')
