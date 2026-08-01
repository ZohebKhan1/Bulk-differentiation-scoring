# ----
# author:
# - Zoheb Khan
#
# script path:
# - tests/test_score_differentiation_timing.R
#
# functions:
# - functions/score_differentiation_timing.R
# ----

source('functions/score_differentiation_timing.R')

expect_error_matching <- function(code, pattern) {
  tryCatch(
    {
      force(code)
      FALSE
    },
    error = function(error_condition) {
      grepl(pattern, conditionMessage(error_condition), fixed = TRUE)
    }
  )
}

expression_matrix <- rbind(
  gene_1 = c(0, 0.2, 1.0, 1.2, 2.0, 2.2),
  gene_2 = c(2.3, 2.1, 1.3, 1.1, 0.3, 0.1),
  gene_3 = c(0.1, 0.2, 1.2, 1.0, 0.4, 0.5),
  gene_4 = c(0.3, 0.4, 0.6, 0.8, 1.4, 1.2)
)
colnames(expression_matrix) <- paste0('sample_', seq_len(ncol(expression_matrix)))

metadata <- data.frame(
  sample_id = colnames(expression_matrix),
  day = rep(1:3, each = 2),
  cohort = c(rep('reference', 4), rep('test', 2)),
  stringsAsFactors = FALSE
)
metadata <- metadata[c(4, 2, 6, 1, 5, 3), ]

all_reference_fit <- score_differentiation_timing(
  expression_matrix = expression_matrix,
  metadata = metadata,
  temporal_genes = rownames(expression_matrix),
  sample_id_col = 'sample_id',
  time_col = 'day',
  n_pcs = 3
)

stopifnot(
  identical(all_reference_fit$scores$sample_id, colnames(expression_matrix)),
  all(all_reference_fit$scores$is_reference),
  all(all_reference_fit$scores$predicted_time >= 1),
  all(all_reference_fit$scores$predicted_time <= 3),
  all(all_reference_fit$scores$differentiation_score >= 0),
  all(all_reference_fit$scores$differentiation_score <= 1),
  all(is.finite(all_reference_fit$scores$squared_distance))
)

heldout_fit <- score_differentiation_timing(
  expression_matrix = expression_matrix,
  metadata = metadata,
  temporal_genes = rownames(expression_matrix),
  sample_id_col = 'sample_id',
  time_col = 'day',
  reference_col = 'cohort',
  reference_values = 'reference',
  n_pcs = 2
)

stopifnot(
  identical(heldout_fit$scores$is_reference, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)),
  all(heldout_fit$scores$predicted_time >= 1),
  all(heldout_fit$scores$predicted_time <= 2),
  all(heldout_fit$scores$differentiation_score >= 0),
  all(heldout_fit$scores$differentiation_score <= 1)
)

invalid_metadata <- rbind(metadata, metadata[1, , drop = FALSE])
stopifnot(
  expect_error_matching(
    score_differentiation_timing(
      expression_matrix = expression_matrix,
      metadata = invalid_metadata,
      temporal_genes = rownames(expression_matrix),
      time_col = 'day'
    ),
    'unique'
  )
)

nonfinite_metadata <- metadata
nonfinite_metadata$day[nonfinite_metadata$cohort == 'test'] <- Inf
stopifnot(
  expect_error_matching(
    score_differentiation_timing(
      expression_matrix = expression_matrix,
      metadata = nonfinite_metadata,
      temporal_genes = rownames(expression_matrix),
      time_col = 'day'
    ),
    'finite numeric'
  )
)

message('Standalone scorer contract checks passed.')
