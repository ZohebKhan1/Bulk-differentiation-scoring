library(DESeq2)
library(edgeR)

# source project functions
source('functions/load_GSE122380_data.R')
source('functions/select_temporal_genes.R')
source('functions/score_differentiation_timing.R')

# 1.0 define script parameters and paths -----------------

expression_cpm_cutoff = 10
lrt_padj_cutoff = 1e-7
vst_dynamic_range_cutoff = 0.6
n_pcs = 3L

validation_results_path = 'tmp/GSE122380_leave_one_cell_line_out_validation.rds'

# 1.1 read direct inputs -----------------

GSE122380_data <- load_GSE122380_data()
metadata <- GSE122380_data$metadata
counts <- GSE122380_data$counts
vst <- GSE122380_data$vst

ordered_days <- sort(unique(metadata$day_numeric))
heldout_cell_lines <- sort(unique(as.character(metadata$cell_line)))

# 1.2 define local helpers -----------------

split_temporal_genes <- function(temporal_genes, training_metadata) {
  gene_day_means <- sapply(ordered_days, function(day_value) {
    day_samples <- training_metadata$sample_id[
      training_metadata$day_numeric == day_value
    ]
    rowMeans(vst[temporal_genes, day_samples, drop = FALSE], na.rm = TRUE)
  })
  gene_day_correlation <- apply(gene_day_means, 1, function(gene_values) {
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
  scoring_metadata <- metadata
  scoring_metadata$is_training <- scoring_metadata$sample_id %in%
    training_metadata$sample_id

  timing_fit <- score_differentiation_timing(
    expression_matrix = vst,
    metadata = scoring_metadata,
    temporal_genes = gene_ids,
    sample_id_col = 'sample_id',
    time_col = 'day_numeric',
    reference_col = 'is_training',
    reference_values = TRUE,
    n_pcs = n_pcs
  )
  heldout_scores <- timing_fit$scores[!timing_fit$scores$is_reference, , drop = FALSE]
  pca_variance <- summary(timing_fit$pca_fit)$importance[2, seq_len(n_pcs)] * 100

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
    gene_count = length(gene_ids),
    pca_variance_percent = round(pca_variance, 2),
    scores = score_table
  )
}

analyze_heldout_cell_line <- function(heldout_cell_line) {
  message('Analyzing held-out cell line: ', heldout_cell_line)
  training_metadata <- metadata[
    metadata$cell_line != heldout_cell_line,
    ,
    drop = FALSE
  ]
  training_metadata$cell_line <- droplevels(training_metadata$cell_line)
  training_sample_ids <- training_metadata$sample_id

  temporal_selection <- select_temporal_genes(
    counts = counts[, training_sample_ids, drop = FALSE],
    vst = vst[, training_sample_ids, drop = FALSE],
    metadata = training_metadata,
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
    data.frame(
      gene_set = result$gene_set,
      gene_count = result$gene_count,
      pc1_percent = result$pca_variance_percent[[1]],
      pc2_percent = result$pca_variance_percent[[2]],
      pc3_percent = result$pca_variance_percent[[3]],
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

# 2.0 run cell-line cross-validation -----------------

validation_results <- setNames(
  lapply(heldout_cell_lines, analyze_heldout_cell_line),
  heldout_cell_lines
)
score_table <- do.call(rbind, lapply(validation_results, function(result) result$scores))
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
      'pc1_percent',
      'pc2_percent',
      'pc3_percent'
    ),
    drop = FALSE
  ]
}))
rownames(summary_table) <- NULL
validation_output <- list(
  scores = score_table,
  summary = summary_table,
  results = validation_results
)
dir.create(dirname(validation_results_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(validation_output, validation_results_path)

# 3.0 report summary -----------------

validation_output$summary
validation_results_path
