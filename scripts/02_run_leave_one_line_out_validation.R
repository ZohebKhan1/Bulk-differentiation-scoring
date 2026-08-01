# ----
# author:
# - Zoheb Khan
#
# script path:
# - scripts/02_run_leave_one_line_out_validation.R
#
# functions:
# - functions/load_GSE122380_data.R
# - functions/select_temporal_genes.R
# - functions/score_differentiation_timing.R
#
# params:
# - config/analysis.yml
#
# input data:
# - data/GSE122380_metadata.rds
# - data/GSE122380_counts.rds
# - data/GSE122380_vst.rds
#
# outputs:
# - tmp/GSE122380_leave_one_line_out_validation.rds
# ----

# 0.0 source packages and dependencies -----------------

required_packages <- c('DESeq2', 'edgeR', 'yaml')
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    'Missing required packages: ',
    paste(missing_packages, collapse = ', '),
    '.',
    call. = FALSE
  )
}

# source canonical loader, temporal selection, and timing scorer
source('functions/load_GSE122380_data.R')
source('functions/select_temporal_genes.R')
source('functions/score_differentiation_timing.R')

# 1.0 read shared parameters -----------------

analysis_params <- yaml::read_yaml('config/analysis.yml')

# 1.1 define script parameters and paths -----------------

expression_cpm_cutoff <- analysis_params$temporal_gene_selection$expression_cpm_cutoff
lrt_padj_cutoff <- analysis_params$temporal_gene_selection$lrt_padj_cutoff
vst_dynamic_range_cutoff <- analysis_params$temporal_gene_selection$vst_dynamic_range_cutoff
n_pcs <- analysis_params$timing_score$n_pcs

validation_cache_path <- 'tmp/GSE122380_leave_one_line_out_validation.rds'

# 1.2 read direct inputs -----------------

GSE122380_data <- load_GSE122380_data()
metadata <- GSE122380_data$metadata
counts <- GSE122380_data$counts
vst <- GSE122380_data$vst

ordered_days <- sort(unique(metadata$day_numeric))
selected_cell_lines <- sort(unique(as.character(metadata$cell_line)))

cache_source_paths <- c(
  'scripts/02_run_leave_one_line_out_validation.R',
  'functions/load_GSE122380_data.R',
  'functions/select_temporal_genes.R',
  'functions/score_differentiation_timing.R',
  'config/analysis.yml'
)
validation_cache_key <- list(
  input_md5 = GSE122380_data$input_md5,
  implementation_md5 = unname(tools::md5sum(cache_source_paths)),
  package_versions = c(
    DESeq2 = as.character(utils::packageVersion('DESeq2')),
    edgeR = as.character(utils::packageVersion('edgeR'))
  ),
  expression_cpm_cutoff = expression_cpm_cutoff,
  lrt_padj_cutoff = lrt_padj_cutoff,
  vst_dynamic_range_cutoff = vst_dynamic_range_cutoff,
  n_pcs = n_pcs,
  heldout_lines = selected_cell_lines
)

# 1.3 define local helpers -----------------

read_validation_cache <- function() {
  if (!file.exists(validation_cache_path)) {
    return(NULL)
  }

  cache <- readRDS(validation_cache_path)
  if (!isTRUE(identical(cache$cache_key, validation_cache_key))) {
    return(NULL)
  }
  cache
}

write_validation_cache <- function(results) {
  score_table <- do.call(rbind, lapply(results, function(result) result$scores))
  rownames(score_table) <- NULL

  summary_table <- do.call(rbind, lapply(results, function(result) {
    summary <- result$pca_summary
    summary$heldout_line <- result$heldout_line
    summary$temporal_genes <- length(result$temporal_genes)
    summary[
      ,
      c(
        'heldout_line',
        'gene_set',
        'gene_count',
        'temporal_genes',
        'pc1_percent',
        'pc2_percent',
        'pc3_percent'
      ),
      drop = FALSE
    ]
  }))
  rownames(summary_table) <- NULL

  validation_cache <- list(
    scores = score_table,
    summary = summary_table,
    results = results,
    selected_cell_lines = names(results),
    cache_key = validation_cache_key,
    notes = list(
      validation_scope = paste(
        'Cell-line cross-validation refits temporal selection, PCA, and the',
        'polyline within the processed GSE122380 cohort; the tracked VST',
        'matrix is shared across folds.'
      ),
      gene_set_split = paste(
        'Maturation and progenitor genes have positive and negative',
        'training-set Spearman correlations with day, respectively.'
      )
    )
  )

  dir.create(dirname(validation_cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(validation_cache, validation_cache_path)
  validation_cache
}

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
  heldout_line,
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
  heldout_scores <- heldout_scores[
    heldout_scores$sample_id %in% metadata$sample_id[metadata$cell_line == heldout_line],
    ,
    drop = FALSE
  ]
  pca_variance <- summary(timing_fit$pca_fit)$importance[2, seq_len(n_pcs)] * 100

  score_table <- data.frame(
    heldout_line = heldout_line,
    sample_id = heldout_scores$sample_id,
    cell_line = heldout_line,
    actual_day = heldout_scores$observed_time,
    gene_set = gene_set_name,
    method = 'Polyline',
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

analyze_heldout_line <- function(heldout_line) {
  message('Analyzing held-out line: ', heldout_line)
  training_metadata <- metadata[metadata$cell_line != heldout_line, , drop = FALSE]
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
  if (any(vapply(gene_sets, length, integer(1)) < 2L)) {
    stop(
      'Held-out line ', heldout_line,
      ' produced a temporal gene subset with fewer than two genes.',
      call. = FALSE
    )
  }

  gene_set_results <- Map(
    score_gene_set,
    gene_sets,
    names(gene_sets),
    MoreArgs = list(
      heldout_line = heldout_line,
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
    heldout_line = heldout_line,
    temporal_genes = temporal_selection$temporal_genes,
    temporal_selection_summary = temporal_selection$summary,
    gene_sets = gene_sets,
    pca_summary = pca_summary,
    scores = score_table
  )
}

# 2.0 run cell-line cross-validation -----------------

validation_cache <- read_validation_cache()
if (is.null(validation_cache)) {
  validation_results <- list()
} else {
  validation_results <- validation_cache$results
}

for (heldout_line in selected_cell_lines) {
  if (!is.null(validation_results[[heldout_line]])) {
    message('Using cached result for held-out line: ', heldout_line)
    next
  }

  validation_results[[heldout_line]] <- analyze_heldout_line(heldout_line)
  validation_results <- validation_results[
    selected_cell_lines[selected_cell_lines %in% names(validation_results)]
  ]
  validation_cache <- write_validation_cache(validation_results)
}

validation_results <- validation_results[selected_cell_lines]
validation_cache <- write_validation_cache(validation_results)

# 3.0 report summary -----------------

validation_cache$summary
validation_cache_path
