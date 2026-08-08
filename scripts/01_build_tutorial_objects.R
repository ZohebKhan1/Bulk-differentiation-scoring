library(AnnotationDbi)
library(circlize)
library(clusterProfiler)
library(ComplexHeatmap)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(patchwork)
library(ragg)
library(scales)
library(systemfonts)
library(viridis)

# source project functions
source('src/load_GSE122380_data.R', local = TRUE)
source('src/run_leave_one_cell_line_out_validation.R', local = TRUE)
source('functions/select_temporal_genes.R', local = TRUE)
source('functions/score_differentiation_timing.R', local = TRUE)

# 1.0 define script parameters and paths -----------------

expression_cpm_cutoff = 10
lrt_padj_cutoff = 1e-7
vst_dynamic_range_cutoff = 0.6

docs_figure_dir = 'docs/assets/figures'
tutorial_font_dir = 'tutorial/assets/fonts'
loo_validation_cache_path = 'cache/GSE122380_leave_one_cell_line_out_validation.rds'

reference_overview_figure_width = 7.20
reference_overview_figure_height = 3.81

temporal_heatmap_figure_width = 7.20
temporal_heatmap_figure_height = 4.42

temporal_clusters_figure_width = 7.20
temporal_clusters_figure_height = 8.52

pca_day_figure_width = 7.20
pca_day_figure_height = 8.52

pc1_validation_figure_width = 7.20
pc1_validation_figure_height = 4.95

timing_polyline_figure_width = 7.20
timing_polyline_figure_height = 3.81

score_by_day_figure_width = 7.20
score_by_day_figure_height = 3.65

loo_cell_line_predictions_figure_width = 7.20
loo_cell_line_predictions_figure_height = 5.75

loo_summary_figure_width = 7.20
loo_summary_figure_height = 3.15

loo_timepoint_accuracy_figure_width = 6.20
loo_timepoint_accuracy_figure_height = 3.10

figure_dpi = 600
figure_family = 'GSE122380 Nimbus Sans'
# ggplot2 converts text from 72.27 typographic points per inch to millimetres
ggplot_points_per_mm = 72.27 / 25.4

n_heatmap_genes = 1500L
n_temporal_clusters = 4L
n_pc1_loading_plot_genes_per_direction = 10L
n_pc1_go_genes_per_direction = 500L
pca_gene_fraction = 0.10
pc1_negative_color = '#1E40AF'
pc1_positive_color = '#A80000'
go_top_n_terms = 5L
go_term_padj_cutoff = 0.05
go_min_gene_count_in_term = 10L
go_min_genes_in_go_db = 26L
go_max_genes_in_go_db = 499L
reference_day_palette <- viridis::viridis(15, option = 'D')
reference_day_palette[[15]] <- '#D8B11E'
annotation_day_palette <- grDevices::colorRampPalette(reference_day_palette)(256)
correlation_palette <- c(
  '#093F60', '#176086', '#2C83AA', '#56A5B8', '#82B6BB',
  '#AECFC0', '#D5E3BB', '#F6E699', '#FAD171', '#F5B14A',
  '#EA832A', '#D95F24', '#C43C22', '#A92325', '#831026'
)

systemfonts::register_font(
  name = figure_family,
  plain = file.path(tutorial_font_dir, 'NimbusSans-Regular.otf'),
  bold = file.path(tutorial_font_dir, 'NimbusSans-Bold.otf'),
  italic = file.path(tutorial_font_dir, 'NimbusSans-Italic.otf'),
  bolditalic = file.path(tutorial_font_dir, 'NimbusSans-BoldItalic.otf')
)
# ComplexHeatmap measures legends on an internal PDF device. Nimbus Sans is
# Helvetica-compatible, so the same fallback metrics preserve its layout.
base::do.call(
  grDevices::pdfFonts,
  stats::setNames(grDevices::pdfFonts('Helvetica'), figure_family)
)

# 1.2 read direct inputs -----------------

GSE122380_data <- load_GSE122380_data()
metadata <- GSE122380_data$metadata
counts <- GSE122380_data$counts
vst <- GSE122380_data$vst

# 1.3 define local helpers -----------------

text_size <- function(size_in_points) {
  size_in_points / ggplot_points_per_mm
}

theme_tutorial <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = 'black', family = figure_family),
      axis.title = ggplot2::element_text(color = 'black', family = figure_family, size = 7),
      axis.text = ggplot2::element_text(color = 'black', family = figure_family, size = 6),
      axis.line = ggplot2::element_line(color = 'black', linewidth = 0.24),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = 0.24),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'bold',
        size = 6
      ),
      legend.text = ggplot2::element_text(color = 'black', family = figure_family, size = 5.5),
      strip.text = ggplot2::element_text(color = 'black', family = figure_family, size = 6),
      plot.title = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'plain',
        size = 7
      ),
      plot.subtitle = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'plain',
        size = 6
      ),
      plot.tag = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'bold',
        size = 8
      )
    )
}

scale_color_by_day <- function(name = 'Day', option = 'D') {
  ggplot2::scale_color_viridis_c(
    name = name,
    option = option,
    breaks = c(1, 5, 10, 15),
    guide = ggplot2::guide_colorbar(
      direction = 'horizontal',
      title.position = 'top',
      title.hjust = 0.5,
      barwidth = grid::unit(0.78, 'in'),
      barheight = grid::unit(0.06, 'in'),
      frame.colour = 'black',
      frame.linewidth = 0.18,
      ticks = FALSE,
      theme = ggplot2::theme(
        legend.title = ggplot2::element_text(size = 6, face = 'bold', family = figure_family),
        legend.text = ggplot2::element_text(size = 5.5, family = figure_family),
        legend.ticks = ggplot2::element_blank(),
        legend.ticks.length = grid::unit(0, 'pt')
      )
    )
  )
}

theme_legend_overlay <- function() {
  ggplot2::theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = 'horizontal',
    legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key.width = grid::unit(0.30, 'in')
  )
}

theme_legend_top <- function() {
  ggplot2::theme(
    legend.position = 'top',
    legend.justification = 'center',
    legend.direction = 'horizontal',
    legend.background = ggplot2::element_blank(),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 1, 0),
    legend.key.width = grid::unit(0.30, 'in')
  )
}

save_figure <- function(path_stub, plot, width, height) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = paste0(path_stub, '.png'),
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    dpi = figure_dpi,
    bg = 'white',
    limitsize = FALSE
  )
}

save_drawn_figure <- function(path_stub, draw_fn, width, height, png_scale = 1) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(
    filename = paste0(path_stub, '.png'),
    width = width * png_scale,
    height = height * png_scale,
    units = 'in',
    res = figure_dpi,
    background = 'white'
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  draw_fn()
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

label_ensembl_genes <- function(gene_ids) {
  clean_ids <- sub('\\..*$', '', gene_ids)
  symbol_lists <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = clean_ids,
    keytype = 'ENSEMBL',
    column = 'SYMBOL',
    multiVals = 'CharacterList'
  )

  vapply(seq_along(gene_ids), function(gene_index) {
    symbols <- unique(as.character(symbol_lists[[clean_ids[[gene_index]]]]))
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    if (length(symbols) == 1L) symbols else gene_ids[[gene_index]]
  }, character(1))
}

orient_pca_fit_by_day <- function(pca_fit, day_by_sample) {
  pc1_day_cor <- stats::cor(
    pca_fit$x[, 'PC1'],
    day_by_sample[match(rownames(pca_fit$x), names(day_by_sample))],
    use = 'pairwise.complete.obs'
  )
  if (is.finite(pc1_day_cor) && pc1_day_cor < 0) {
    pca_fit$x[, 'PC1'] <- -pca_fit$x[, 'PC1']
    pca_fit$rotation[, 'PC1'] <- -pca_fit$rotation[, 'PC1']
  }
  pca_fit
}

get_pca_color_option <- function(set_label) {
  switch(
    set_label,
    'C1+C2 (Early)' = 'inferno',
    'C3+C4 (Late)' = 'mako',
    'viridis'
  )
}

format_correlation_label <- function(x, y) {
  r_value <- stats::cor(x, y, use = 'pairwise.complete.obs')
  sprintf('r = %.2f\nr\u00b2 = %.2f', r_value, r_value^2)
}

format_go_label <- function(x, width = 34L) {
  vapply(
    x,
    function(term) paste(strwrap(term, width = width), collapse = '\n'),
    character(1)
  )
}

make_go_count_legend_breaks <- function(counts) {
  as.numeric(stats::quantile(
    counts,
    probs = c(0, 0.25, 0.75, 1),
    type = 7,
    names = FALSE
  ))
}

format_go_count_labels <- function(counts) {
  ifelse(
    abs(counts - round(counts)) < 0.01,
    as.character(round(counts)),
    formatC(counts, format = 'f', digits = 1)
  )
}

run_go_enrichment <- function(gene_ids, universe_gene_ids) {
  clean_gene_ids <- unique(sub('\\..*$', '', gene_ids))
  clean_universe <- unique(sub('\\..*$', '', universe_gene_ids))

  go_fit <- clusterProfiler::enrichGO(
    gene = clean_gene_ids,
    universe = clean_universe,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = go_term_padj_cutoff,
    qvalueCutoff = 1,
    minGSSize = go_min_genes_in_go_db,
    maxGSSize = go_max_genes_in_go_db,
    readable = TRUE
  )

  go_results <- as.data.frame(go_fit)
  if (nrow(go_results) == 0L) {
    return(go_results)
  }

  go_results <- go_results[
    !is.na(go_results$p.adjust) &
      go_results$p.adjust < go_term_padj_cutoff &
      go_results$Count >= go_min_gene_count_in_term,
    ,
    drop = FALSE
  ]
  if (nrow(go_results) == 0L) {
    return(go_results)
  }

  go_results$neg_log10_pvalue <- -log10(go_results$pvalue)
  go_results <- go_results[
    order(go_results$neg_log10_pvalue, go_results$Count, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  utils::head(go_results, go_top_n_terms)
}

make_go_dotplot <- function(go_results, plot_title, point_color) {
  if (nrow(go_results) == 0L) {
    empty_df <- data.frame(x = 0, y = plot_title)
    return(
      ggplot2::ggplot(empty_df, ggplot2::aes(x, y)) +
        ggplot2::geom_blank() +
        ggplot2::annotate(
          'text',
          x = 0,
          y = plot_title,
          label = 'No enriched GO terms',
          family = figure_family,
          size = text_size(5.5)
        ) +
        ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
        theme_tutorial() +
        ggplot2::theme(
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(
            color = point_color,
            face = 'bold',
            family = figure_family,
            size = 7
          )
        )
    )
  }

  plot_df <- go_results
  plot_df$description_label <- format_go_label(plot_df$Description)
  plot_df$description_label <- factor(
    plot_df$description_label,
    levels = rev(plot_df$description_label)
  )
  count_legend_breaks <- make_go_count_legend_breaks(plot_df$Count)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(neg_log10_pvalue, description_label)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = Count),
      color = point_color,
      alpha = 0.96
    ) +
    ggplot2::scale_size_continuous(
      name = '# genes',
      breaks = count_legend_breaks,
      labels = format_go_count_labels(count_legend_breaks),
      range = c(1.1, 3.0),
      guide = ggplot2::guide_legend(
        direction = 'vertical',
        title.position = 'top',
        title.hjust = 0,
        override.aes = list(alpha = 0.96)
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      x = '-log10(p)',
      y = NULL
    ) +
    theme_tutorial() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = point_color,
        face = 'bold',
        family = figure_family,
        size = 7
      ),
      axis.text.y = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        lineheight = 0.9,
        size = 5.5
      ),
      legend.position = 'right',
      legend.direction = 'vertical',
      legend.title = ggplot2::element_text(face = 'bold', family = figure_family, size = 6),
      legend.text = ggplot2::element_text(family = figure_family, size = 5.5),
      legend.key.size = grid::unit(0.26, 'cm'),
      legend.margin = ggplot2::margin(0, 0, 0, 2),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(4, 3, 4, 3)
    )
}

format_regression_label <- function(x, y) {
  lm_fit <- stats::lm(y ~ x)
  coefs <- stats::coef(lm_fit)
  r_value <- stats::cor(x, y, use = 'pairwise.complete.obs')
  intercept <- coefs[[1]]
  slope <- coefs[[2]]
  sign_text <- if (slope < 0) '-' else '+'
  sprintf(
    'r = %.2f\nr\u00b2 = %.2f\ny = %.2f %s %.2fx',
    r_value,
    r_value^2,
    intercept,
    sign_text,
    abs(slope)
  )
}

plot_reference_pca <- function(pca_data, day_path, x_label, y_label, pca_gene_count) {
  reference_day_colors <- stats::setNames(reference_day_palette, seq_len(15))
  phase_colors <- c(
    'Pluripotent' = reference_day_colors[['1']],
    'Mesoderm' = reference_day_colors[['4']],
    'Immature cardiomyocyte' = reference_day_colors[['8']],
    'Mature cardiomyocyte' = reference_day_colors[['15']]
  )
  phase_label_fills <- phase_colors
  phase_label_fills[['Mature cardiomyocyte']] <- '#B8860B'

  pca_data$phase <- ifelse(
    pca_data$day_numeric <= 2,
    'Pluripotent',
    ifelse(
      pca_data$day_numeric <= 5,
      'Mesoderm',
      ifelse(
        pca_data$day_numeric <= 9,
        'Immature cardiomyocyte',
        'Mature cardiomyocyte'
      )
    )
  )
  pca_data$phase <- factor(pca_data$phase, levels = names(phase_colors))

  ellipse_df <- base::do.call(rbind, lapply(names(phase_colors), function(phase_name) {
    phase_data <- pca_data[pca_data$phase == phase_name, c('PC1', 'PC2'), drop = FALSE]
    center <- colMeans(phase_data)
    covariance <- stats::cov(phase_data)
    ellipse_angle <- seq(0, 2 * pi, length.out = 160)
    ellipse_circle <- rbind(cos(ellipse_angle), sin(ellipse_angle))
    ellipse_coords <- t(
      center +
        sqrt(stats::qchisq(0.80, df = 2)) *
          t(chol(covariance)) %*% ellipse_circle
    )

    data.frame(
      phase = phase_name,
      PC1 = ellipse_coords[, 1],
      PC2 = ellipse_coords[, 2],
      stringsAsFactors = FALSE
    )
  }))

  ellipse_layers <- lapply(names(phase_colors), function(phase_name) {
    ggplot2::geom_path(
      data = ellipse_df[ellipse_df$phase == phase_name, , drop = FALSE],
      ggplot2::aes(PC1, PC2),
      inherit.aes = FALSE,
      color = phase_colors[[phase_name]],
      linewidth = 0.50,
      alpha = 0.9
    )
  })

  day_label_offsets <- data.frame(
    day_numeric = seq_len(15),
    nudge_x = c(-2, -4, -2, 0, 0, 0, -2, -2, -8, 5, -10, 8, -9, 10, 0),
    nudge_y = c(-6, 5, -5, 5, -5, 6, -7, 6, -5, 5, -9, 0, -2, -5, -12)
  )
  day_path <- merge(day_path, day_label_offsets, by = 'day_numeric', all.x = TRUE)
  day_path <- day_path[order(day_path$day_numeric), ]
  day_path$day_label <- paste0('D', day_path$day_numeric)
  day_path$label_x <- day_path$PC1 + day_path$nudge_x
  day_path$label_y <- day_path$PC2 + day_path$nudge_y

  phase_labels <- data.frame(
    day_numeric = c(1, 4, 8, 14),
    phase = c('Pluripotent', 'Mesoderm', 'Immature cardiomyocyte', 'Mature cardiomyocyte'),
    nudge_x = c(-4, -22, 8, 2),
    nudge_y = c(-12, 4, 19, -16),
    stringsAsFactors = FALSE
  )
  phase_labels <- merge(
    phase_labels,
    day_path[, c('day_numeric', 'PC1', 'PC2')],
    by = 'day_numeric',
    all.x = TRUE
  )
  phase_labels$label_x <- phase_labels$PC1 + phase_labels$nudge_x
  phase_labels$label_y <- phase_labels$PC2 + phase_labels$nudge_y
  phase_labels$phase <- factor(phase_labels$phase, levels = names(phase_colors))

  plot_bounds <- rbind(
    pca_data[, c('PC1', 'PC2')],
    ellipse_df[, c('PC1', 'PC2')],
    stats::setNames(day_path[, c('label_x', 'label_y')], c('PC1', 'PC2')),
    stats::setNames(phase_labels[, c('label_x', 'label_y')], c('PC1', 'PC2'))
  )
  pca_xlim <- range(plot_bounds$PC1, na.rm = TRUE) + c(-8, 8)
  pca_ylim <- range(plot_bounds$PC2, na.rm = TRUE) + c(-8, 8)

  ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2)) +
    ellipse_layers +
    ggplot2::geom_point(
      ggplot2::aes(color = day_numeric),
      size = 0.95,
      alpha = 0.92
    ) +
    ggplot2::geom_text(
      data = day_path,
      ggplot2::aes(label_x, label_y, label = day_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = text_size(5.5),
      color = 'black'
    ) +
    ggplot2::geom_label(
      data = phase_labels,
      ggplot2::aes(label_x, label_y, label = phase, fill = phase),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = text_size(5.5),
      color = 'white',
      linewidth = 0,
      label.r = grid::unit(0.08, 'lines'),
      label.padding = grid::unit(0.18, 'lines')
    ) +
    ggplot2::scale_color_gradientn(
      colors = reference_day_palette,
      name = 'Differentiation day',
      breaks = c(1, 5, 10, 15),
      guide = ggplot2::guide_colorbar(
        title.position = 'top',
        title.hjust = 0,
        barwidth = grid::unit(0.78, 'in'),
        barheight = grid::unit(0.06, 'in'),
        frame.colour = 'black',
        frame.linewidth = 0.18,
        theme = ggplot2::theme(
          legend.title = ggplot2::element_text(family = figure_family),
          legend.text = ggplot2::element_text(family = figure_family),
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')
        )
      )
    ) +
    ggplot2::scale_fill_manual(values = phase_label_fills, guide = 'none') +
    ggplot2::labs(
      title = 'PCA: GSE122380',
      subtitle = paste0('Top 10% variable genes (n=', scales::comma(pca_gene_count), ')'),
      tag = 'a',
      x = x_label,
      y = y_label
    ) +
    ggplot2::coord_cartesian(xlim = pca_xlim, ylim = pca_ylim, clip = 'off') +
    ggplot2::theme_classic(base_size = 7, base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      plot.title = ggplot2::element_text(family = figure_family, face = 'plain', size = 7),
      plot.subtitle = ggplot2::element_text(family = figure_family, face = 'plain', size = 6),
      plot.tag = ggplot2::element_text(family = figure_family, face = 'bold', size = 8),
      plot.tag.position = c(0.005, 0.995),
      axis.title = ggplot2::element_text(family = figure_family, size = 7),
      axis.text = ggplot2::element_text(family = figure_family, size = 6),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = 'horizontal',
      legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
      legend.margin = ggplot2::margin(1, 1, 1, 1),
      legend.title = ggplot2::element_text(family = figure_family, size = 6, face = 'bold'),
      legend.text = ggplot2::element_text(family = figure_family, size = 5.5),
      legend.key.width = grid::unit(0.30, 'in'),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      plot.margin = ggplot2::margin(10, 8, 6, 10)
    )
}

make_timepoint_correlation_matrix <- function(expression_mat, metadata) {
  day_order <- sort(unique(metadata$day_numeric))
  day_mean_mat <- sapply(day_order, function(day_value) {
    sample_ids <- metadata$sample_id[metadata$day_numeric == day_value]
    rowMeans(expression_mat[, sample_ids, drop = FALSE], na.rm = TRUE)
  })
  colnames(day_mean_mat) <- paste0('D', day_order)

  stats::cor(day_mean_mat, method = 'pearson', use = 'pairwise.complete.obs')
}

get_day_colors <- function(days_to_color) {
  day_index <- round(scales::rescale(
    days_to_color,
    to = c(1, 256),
    from = range(days_to_color)
  ))
  stats::setNames(annotation_day_palette[day_index], paste0('D', days_to_color))
}

plot_timepoint_correlation_heatmap <- function(correlation_mat) {
  ordered_labels <- paste0('D', seq_len(ncol(correlation_mat)))
  y_labels <- rev(ordered_labels)
  n_labels <- length(ordered_labels)
  gap_after <- ceiling(n_labels / 2)
  gap_size <- 0.08
  add_heatmap_gap <- function(index) {
    index + ifelse(index > gap_after, gap_size, 0)
  }
  x_positions <- add_heatmap_gap(seq_along(ordered_labels))
  y_positions <- add_heatmap_gap(seq_along(y_labels))
  max_x_position <- max(x_positions)
  max_y_position <- max(y_positions)

  heatmap_data <- as.data.frame(as.table(correlation_mat[ordered_labels, y_labels]))
  names(heatmap_data) <- c('x_label', 'y_label', 'correlation')
  heatmap_data$x_index <- add_heatmap_gap(match(heatmap_data$x_label, ordered_labels))
  heatmap_data$y_index <- add_heatmap_gap(match(heatmap_data$y_label, y_labels))

  ordered_days <- as.integer(sub('^D', '', ordered_labels))
  y_days <- as.integer(sub('^D', '', y_labels))
  day_colors <- get_day_colors(sort(unique(ordered_days)))
  top_annotation <- data.frame(
    x_index = x_positions,
    y_index = max_y_position + 0.82,
    fill = unname(day_colors[ordered_labels])
  )
  left_annotation <- data.frame(
    x_index = 0.15,
    y_index = y_positions,
    fill = unname(day_colors[paste0('D', y_days)])
  )
  top_labels <- data.frame(
    x_index = x_positions,
    y_index = max_y_position + 1.37,
    label = ordered_labels
  )

  ggplot2::ggplot() +
    ggplot2::annotate(
      'segment',
      x = 1,
      xend = max_x_position,
      y = max_y_position + 2.45,
      yend = max_y_position + 2.45,
      linewidth = 0.20,
      arrow = grid::arrow(length = grid::unit(0.055, 'in'), type = 'closed')
    ) +
    ggplot2::annotate(
      'text',
      x = (max_x_position + 1) / 2,
      y = max_y_position + 2.84,
      label = 'Developmental time',
      family = figure_family,
      fontface = 'bold',
      size = text_size(5.5)
    ) +
    ggplot2::geom_tile(
      data = heatmap_data,
      ggplot2::aes(x_index, y_index, fill = correlation),
      width = 1.01,
      height = 1.01,
      color = NA,
      linewidth = 0
    ) +
    ggplot2::geom_tile(
      data = top_annotation,
      ggplot2::aes(x_index, y_index),
      fill = top_annotation$fill,
      width = 1,
      height = 0.55,
      color = NA
    ) +
    ggplot2::geom_tile(
      data = left_annotation,
      ggplot2::aes(x_index, y_index),
      fill = left_annotation$fill,
      width = 0.55,
      height = 1,
      color = NA
    ) +
    ggplot2::geom_text(
      data = top_labels,
      ggplot2::aes(x_index, y_index, label = label),
      inherit.aes = FALSE,
      family = figure_family,
      size = text_size(5.5),
      angle = 45,
      hjust = 0,
      vjust = 0.5,
      color = 'grey25'
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      position = 'top',
      limits = c(-0.16, max_x_position + 0.51),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_positions,
      labels = y_labels,
      limits = c(0.49, max_y_position + 3.09),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_fill_gradientn(
      colors = correlation_palette,
      limits = c(0, 1),
      name = 'Pearson r',
      breaks = c(0, 0.5, 1),
      labels = function(x) {
        ifelse(
          x %% 1 == 0,
          formatC(x, format = 'f', digits = 0),
          formatC(x, format = 'f', digits = 1)
        )
      },
      guide = ggplot2::guide_colorbar(
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = grid::unit(1.25, 'in'),
        barheight = grid::unit(0.08, 'in'),
        frame.colour = 'black',
        frame.linewidth = 0.18,
        theme = ggplot2::theme(
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')
        )
      )
    ) +
    ggplot2::labs(title = NULL, tag = 'b', x = NULL, y = NULL) +
    ggplot2::coord_fixed(clip = 'off') +
    ggplot2::theme_minimal(base_size = 7, base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(family = figure_family, face = 'bold', size = 8),
      plot.tag.position = c(0.005, 0.995),
      axis.text.x = ggplot2::element_text(
        family = figure_family,
        size = 6,
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        margin = ggplot2::margin(b = -2)
      ),
      axis.text.y = ggplot2::element_text(family = figure_family, size = 6),
      legend.position = 'bottom',
      legend.justification = 'center',
      legend.title = ggplot2::element_text(size = 6, face = 'bold', family = figure_family),
      legend.text = ggplot2::element_text(size = 5.5, family = figure_family),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-24, 0, 0, 34),
      legend.key.width = grid::unit(0.30, 'in'),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      plot.margin = ggplot2::margin(10, 8, 0, 10)
    )
}

make_cluster_trajectory_plot <- function(cluster_name) {
  cluster_df <- cluster_day_gene[cluster_day_gene$cluster == cluster_name, ]
  mean_df <- cluster_means[cluster_means$cluster == cluster_name, ]
  cluster_title <- c(
    C1 = 'Cluster 1 (Stem cell/pluripotent)',
    C2 = 'Cluster 2 (Mesoderm/pluripotent)',
    C3 = 'Cluster 3 (Late cardiomyocyte maturation)',
    C4 = 'Cluster 4 (Early cardiomyocyte maturation)'
  )[[cluster_name]]
  y_limits <- list(
    C1 = c(-2, 3),
    C2 = c(-2, 2.5),
    C3 = c(-2, 2),
    C4 = c(-3, 2)
  )[[cluster_name]]
  y_breaks <- list(
    C1 = c(-2, 0, 3),
    C2 = c(-2, 2),
    C3 = c(-2, 0, 2),
    C4 = c(-3, 0, 2)
  )[[cluster_name]]

  ggplot2::ggplot(cluster_df, ggplot2::aes(day_numeric, zscore, group = gene_id)) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = 'black',
      linetype = 'dashed',
      linewidth = 0.24
    ) +
    ggplot2::geom_line(color = 'gray72', linewidth = 0.18) +
    ggplot2::geom_smooth(
      data = mean_df,
      ggplot2::aes(day_numeric, zscore),
      inherit.aes = FALSE,
      method = 'loess',
      formula = y ~ x,
      se = FALSE,
      span = 0.75,
      color = cluster_colors[[cluster_name]],
      linewidth = 0.75
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(title = cluster_title, x = 'Differentiation day', y = 'Mean VST z-score') +
    theme_tutorial() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = cluster_colors[[cluster_name]],
        face = 'bold',
        family = figure_family,
        size = 7
      )
    )
}

run_cluster_go_enrichment <- function(cluster_assignments, universe_gene_ids) {
  cluster_go <- lapply(levels(cluster_assignments), function(cluster_name) {
    run_go_enrichment(
      gene_ids = names(cluster_assignments)[cluster_assignments == cluster_name],
      universe_gene_ids = universe_gene_ids
    )
  })
  names(cluster_go) <- levels(cluster_assignments)
  cluster_go
}

make_cluster_go_plot <- function(cluster_name) {
  make_go_dotplot(
    go_results = cluster_go_results[[cluster_name]],
    plot_title = paste0(sub('^C', 'Cluster ', cluster_name), ' GO:BP'),
    point_color = cluster_colors[[cluster_name]]
  )
}

calculate_pca <- function(gene_ids, set_label) {
  set_pca_input <- t(vst[gene_ids, metadata$sample_id, drop = FALSE])
  set_day_means <- do.call(rbind, lapply(days, function(day_value) {
    colMeans(set_pca_input[metadata$day_numeric == day_value, , drop = FALSE])
  }))
  set_pca_fit <- stats::prcomp(set_day_means, center = TRUE, scale. = FALSE)
  set_pca_var <- round(summary(set_pca_fit)$importance[2, 1:2] * 100, 2)
  set_pca_coordinates <- sweep(
    set_pca_input,
    MARGIN = 2L,
    STATS = set_pca_fit$center,
    FUN = '-'
  ) %*% set_pca_fit$rotation[, 1:2, drop = FALSE]
  pc1_day_correlation <- stats::cor(set_pca_coordinates[, 'PC1'], metadata$day_numeric)
  if (is.finite(pc1_day_correlation) && pc1_day_correlation < 0) {
    set_pca_coordinates[, 'PC1'] <- -set_pca_coordinates[, 'PC1']
  }

  gene_set_pca_data <- data.frame(
    sample_id = rownames(set_pca_coordinates),
    PC1 = set_pca_coordinates[, 1],
    PC2 = set_pca_coordinates[, 2],
    day_numeric = metadata$day_numeric[
      match(rownames(set_pca_coordinates), metadata$sample_id)
    ],
    gene_set = set_label,
    stringsAsFactors = FALSE
  )

  list(
    data = gene_set_pca_data,
    variance_percent = set_pca_var
  )
}

make_pca_plot <- function(pca_result, plot_title) {
  plot_data <- pca_result$data
  variance_percent <- pca_result$variance_percent
  color_option <- get_pca_color_option(plot_title)

  ggplot2::ggplot(plot_data, ggplot2::aes(PC1, PC2, color = day_numeric)) +
    ggplot2::geom_point(size = 0.95, alpha = 1) +
    scale_color_by_day(name = 'Day', option = color_option) +
    ggplot2::labs(
      title = plot_title,
      x = paste0('PC1 (', variance_percent[[1]], '%)'),
      y = paste0('PC2 (', variance_percent[[2]], '%)')
    ) +
    theme_tutorial() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = figure_family,
        face = 'plain',
        size = 7,
        margin = ggplot2::margin(b = -1)
      ),
      plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5)
    ) +
    theme_legend_overlay()
}

make_pc1_time_plot <- function(pca_result) {
  plot_data <- pca_result$data
  color_option <- get_pca_color_option(plot_data$gene_set[[1]])
  label_text <- format_correlation_label(plot_data$day_numeric, plot_data$PC1)
  label_x <- min(plot_data$day_numeric)
  label_y <- max(plot_data$PC1, na.rm = TRUE)

  ggplot2::ggplot(plot_data, ggplot2::aes(day_numeric, PC1)) +
    ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = 0.95, alpha = 1) +
    ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = 0.45
    ) +
    ggplot2::annotate(
      'text',
      x = label_x,
      y = label_y,
      label = label_text,
      hjust = 0,
      vjust = 1,
      family = figure_family,
      size = text_size(7.5)
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    scale_color_by_day(name = 'Day', option = color_option) +
    ggplot2::labs(
      x = 'Differentiation day',
      y = 'PC1 position'
    ) +
    theme_tutorial() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 14)
    ) +
    theme_legend_top()
}

make_loading_vector_data <- function(component_name) {
  component_loadings <- pca_fit$rotation[, component_name]
  top_pos <- names(
    sort(component_loadings, decreasing = TRUE)
  )[seq_len(n_pc1_loading_plot_genes_per_direction)]
  top_neg <- names(
    sort(component_loadings, decreasing = FALSE)
  )[seq_len(n_pc1_loading_plot_genes_per_direction)]
  top_genes <- c(top_neg, top_pos)

  data.frame(
    gene_id = top_genes,
    gene_label = label_ensembl_genes(top_genes),
    PC1_loading = pca_fit$rotation[top_genes, 'PC1'],
    PC2_loading = pca_fit$rotation[top_genes, 'PC2'],
    direction = rep(
      c(paste0(component_name, '-'), paste0(component_name, '+')),
      each = n_pc1_loading_plot_genes_per_direction
    ),
    stringsAsFactors = FALSE
  )
}

make_loading_vector_plot <- function(loading_df, component_name, plot_title) {
  negative_label <- paste0(component_name, '-')
  positive_label <- paste0(component_name, '+')
  negative_df <- loading_df[loading_df$direction == negative_label, , drop = FALSE]
  positive_df <- loading_df[loading_df$direction == positive_label, , drop = FALSE]
  x_limit <- max(abs(loading_df$PC1_loading)) * 1.45
  y_limit <- max(abs(loading_df$PC2_loading)) * 1.90
  x_limits <- c(-x_limit, x_limit)
  y_limits <- c(-y_limit, y_limit)

  ggplot2::ggplot(loading_df) +
    ggplot2::geom_hline(yintercept = 0, color = 'grey65', linewidth = 0.20) +
    ggplot2::geom_vline(xintercept = 0, color = 'grey65', linewidth = 0.20) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = PC1_loading,
        yend = PC2_loading,
        color = direction
      ),
      linewidth = 0.24,
      arrow = grid::arrow(length = grid::unit(0.045, 'in'), type = 'closed')
    ) +
    ggrepel::geom_text_repel(
      data = negative_df,
      ggplot2::aes(PC1_loading, PC2_loading, label = gene_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'italic',
      size = text_size(5.7),
      color = pc1_negative_color,
      box.padding = grid::unit(0.08, 'lines'),
      point.padding = grid::unit(0.02, 'lines'),
      min.segment.length = 0,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 2026
    ) +
    ggrepel::geom_text_repel(
      data = positive_df,
      ggplot2::aes(PC1_loading, PC2_loading, label = gene_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'italic',
      size = text_size(5.7),
      color = pc1_positive_color,
      box.padding = grid::unit(0.08, 'lines'),
      point.padding = grid::unit(0.02, 'lines'),
      min.segment.length = 0,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 2027
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(
        c(pc1_negative_color, pc1_positive_color),
        c(negative_label, positive_label)
      ),
      guide = 'none'
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'on') +
    ggplot2::labs(
      title = plot_title,
      x = 'PC1 loading',
      y = 'PC2 loading'
    ) +
    theme_tutorial() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = figure_family,
        face = 'bold',
        size = 7,
        hjust = 0.5
      ),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0)),
      plot.margin = ggplot2::margin(4, 3, 4, 3)
    )
}

smooth_score_series <- function(y_values) {
  smooth_days <- seq(min(days), max(days), length.out = 220)
  stats::spline(
    x = score_summary$day_numeric,
    y = y_values,
    xout = smooth_days,
    method = 'natural'
  )$y
}

# 2.0 select temporal genes -----------------

days <- sort(unique(metadata$day_numeric))
day_by_sample <- stats::setNames(metadata$day_numeric, metadata$sample_id)
temporal_selection <- select_temporal_genes(
  counts = counts,
  vst = vst,
  metadata = metadata,
  expression_cpm_cutoff = expression_cpm_cutoff,
  lrt_padj_cutoff = lrt_padj_cutoff,
  vst_dynamic_range_cutoff = vst_dynamic_range_cutoff
)
temporal_genes <- temporal_selection$temporal_genes
lrt_summary <- temporal_selection$summary

# 3.0 create temporal heatmap and clusters -----------------

heatmap_genes <- utils::head(temporal_genes, min(n_heatmap_genes, length(temporal_genes)))
heatmap_expression <- vst[heatmap_genes, metadata$sample_id, drop = FALSE]
heatmap_z <- t(scale(t(heatmap_expression)))
heatmap_z[!is.finite(heatmap_z)] <- 0
heatmap_z[heatmap_z > 2] <- 2
heatmap_z[heatmap_z < -2] <- -2

column_order <- unlist(lapply(days, function(day_value) {
  day_samples <- metadata$sample_id[metadata$day_numeric == day_value]
  day_dist <- stats::dist(t(heatmap_z[, day_samples, drop = FALSE]))
  day_samples[stats::hclust(day_dist, method = 'average')$order]
}), use.names = FALSE)
heatmap_z <- heatmap_z[, column_order, drop = FALSE]
heatmap_metadata <- metadata[match(column_order, metadata$sample_id), , drop = FALSE]

# The four-cluster partition is interpretive and does not affect the primary score
temporal_day_means <- t(vapply(heatmap_genes, function(gene) {
  vapply(days, function(day) {
    samples <- metadata$sample_id[metadata$day_numeric == day]
    mean(vst[gene, samples], na.rm = TRUE)
  }, numeric(1))
}, numeric(length(days))))
colnames(temporal_day_means) <- paste0('D', days)

temporal_day_z <- t(scale(t(temporal_day_means)))
temporal_day_z[!is.finite(temporal_day_z)] <- 0

temporal_smoothing <- lapply(seq_len(nrow(temporal_day_z)), function(gene_index) {
  suppressWarnings({
      loess_fit <- stats::loess(
        temporal_day_z[gene_index, ] ~ days,
        span = 0.55,
        degree = 2,
        family = 'symmetric',
        control = stats::loess.control(surface = 'direct')
      )
      stats::predict(loess_fit, newdata = days)
  })
})
names(temporal_smoothing) <- rownames(temporal_day_z)
temporal_cluster_input <- base::do.call(rbind, temporal_smoothing)
temporal_cluster_input <- t(scale(t(temporal_cluster_input)))
temporal_cluster_tree <- stats::hclust(
  stats::dist(temporal_cluster_input),
  method = 'ward.D2'
)
raw_clusters <- stats::cutree(temporal_cluster_tree, k = n_temporal_clusters)

cluster_day_means <- base::do.call(rbind, lapply(sort(unique(raw_clusters)), function(cl) {
  genes <- names(raw_clusters)[raw_clusters == cl]
  day_values <- colMeans(temporal_day_z[genes, , drop = FALSE])
  data.frame(
    raw_clusters = cl,
    peak_day = days[which.max(day_values)],
    mean_signal = max(day_values),
    stringsAsFactors = FALSE
  )
}))

cluster_order <- cluster_day_means$raw_clusters[
  order(cluster_day_means$peak_day, -cluster_day_means$mean_signal)
]
cluster_map <- stats::setNames(paste0('C', seq_along(cluster_order)), cluster_order)
gene_clusters <- factor(
  cluster_map[as.character(raw_clusters)],
  levels = paste0('C', seq_along(cluster_order))
)
names(gene_clusters) <- names(raw_clusters)

cluster_colors <- stats::setNames(
  c('#0072B2', '#D55E00', '#009E73', '#CC79A7'),
  levels(gene_clusters)
)

heatmap_color_function <- circlize::colorRamp2(
  seq(-2, 2, length.out = 100),
  grDevices::colorRampPalette(correlation_palette)(100)
)

cluster_block_annotation <- ComplexHeatmap::rowAnnotation(
  cluster = ComplexHeatmap::anno_block(
    gp = grid::gpar(fill = unname(cluster_colors[levels(gene_clusters)]), col = NA),
    labels = paste('Cluster', seq_along(levels(gene_clusters))),
    labels_gp = grid::gpar(
      col = 'white',
      fontface = 'bold',
      fontfamily = figure_family,
      fontsize = 12
    ),
    labels_rot = 270,
    width = grid::unit(12, 'mm'),
    show_name = FALSE
  ),
  width = grid::unit(12, 'mm'),
  show_annotation_name = FALSE
)

p_lrt_heatmap <- ComplexHeatmap::Heatmap(
  heatmap_z,
  name = 'Z-score',
  col = heatmap_color_function,
  right_annotation = cluster_block_annotation,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = gene_clusters[rownames(heatmap_z)],
  row_title = NULL,
  row_title_rot = 0,
  row_title_gp = grid::gpar(
    fontsize = 8,
    fontface = 'bold',
    fontfamily = figure_family,
    col = cluster_colors
  ),
  row_gap = grid::unit(1.8, 'mm'),
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = heatmap_metadata$day_factor,
  column_gap = grid::unit(0, 'mm'),
  column_title = paste0('D', levels(heatmap_metadata$day_factor)),
  column_title_gp = grid::gpar(
    fontsize = 15,
    fontface = 'plain',
    fontfamily = figure_family
  ),
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c('-2', '0', '2'),
    tick_length = grid::unit(0, 'mm'),
    border = 'black',
    legend_height = grid::unit(28, 'mm'),
    grid_width = grid::unit(4.8, 'mm'),
    legend_gp = grid::gpar(col = 'black', lwd = 0.5),
    title_gp = grid::gpar(fontface = 'bold', fontfamily = figure_family, fontsize = 10),
    labels_gp = grid::gpar(fontfamily = figure_family, fontsize = 9)
  ),
  border = FALSE,
  rect_gp = grid::gpar(col = NA),
  use_raster = TRUE,
  raster_quality = 24
)

cluster_genes <- names(gene_clusters)
cluster_day_gene <- base::do.call(rbind, lapply(cluster_genes, function(gene) {
  data.frame(
    gene_id = gene,
    cluster = gene_clusters[[gene]],
    day_numeric = days,
    zscore = as.numeric(temporal_day_z[gene, ]),
    stringsAsFactors = FALSE
  )
}))

cluster_means <- stats::aggregate(
  zscore ~ cluster + day_numeric,
  data = cluster_day_gene,
  FUN = mean
)

cluster_go_results <- run_cluster_go_enrichment(
  cluster_assignments = gene_clusters,
  universe_gene_ids = heatmap_genes
)

cluster_panel_plots <- lapply(levels(gene_clusters), function(cluster_name) {
  (make_cluster_trajectory_plot(cluster_name) | make_cluster_go_plot(cluster_name)) +
    patchwork::plot_layout(widths = c(1.9, 0.72))
})

p_cluster_trajectories <- patchwork::wrap_plots(
  cluster_panel_plots,
  ncol = 1
) +
  patchwork::plot_annotation(tag_levels = 'a') &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = 'bold', family = figure_family, size = 8))

# 4.0 create reference dataset overview -----------------

pca_gene_vars <- apply(vst, 1, stats::var, na.rm = TRUE)
pca_gene_count <- ceiling(length(pca_gene_vars) * pca_gene_fraction)
pca_genes <- names(sort(pca_gene_vars, decreasing = TRUE))[seq_len(pca_gene_count)]
pca_input_mat <- vst[pca_genes, metadata$sample_id, drop = FALSE]

overview_pca <- stats::prcomp(t(pca_input_mat), center = TRUE, scale. = FALSE)
overview_pca_df <- data.frame(
  sample_id = rownames(overview_pca$x),
  PC1 = overview_pca$x[, 1],
  PC2 = overview_pca$x[, 2],
  day_numeric = metadata$day_numeric[match(rownames(overview_pca$x), metadata$sample_id)],
  stringsAsFactors = FALSE
)

pc1_flip <- ifelse(
  stats::median(
    overview_pca_df$PC1[overview_pca_df$day_numeric == max(overview_pca_df$day_numeric)],
    na.rm = TRUE
  ) <
    stats::median(
      overview_pca_df$PC1[overview_pca_df$day_numeric == min(overview_pca_df$day_numeric)],
      na.rm = TRUE
    ),
  -1,
  1
)
overview_pca_df$PC1 <- overview_pca_df$PC1 * pc1_flip

overview_day_path <- stats::aggregate(
  overview_pca_df[, c('PC1', 'PC2')],
  by = list(day_numeric = overview_pca_df$day_numeric),
  FUN = stats::median,
  na.rm = TRUE
)
overview_day_path <- overview_day_path[order(overview_day_path$day_numeric), ]
overview_var_explained <- (overview_pca$sdev^2) / sum(overview_pca$sdev^2)
overview_x_label <- sprintf('PC1 (%.1f%%)', overview_var_explained[[1]] * 100)
overview_y_label <- sprintf('PC2 (%.1f%%)', overview_var_explained[[2]] * 100)

p_reference_pca <- plot_reference_pca(
  pca_data = overview_pca_df,
  day_path = overview_day_path,
  x_label = overview_x_label,
  y_label = overview_y_label,
  pca_gene_count = pca_gene_count
)

timepoint_correlation_mat <- make_timepoint_correlation_matrix(
  expression_mat = pca_input_mat,
  metadata = metadata
)
p_reference_timepoint_correlation <- plot_timepoint_correlation_heatmap(
  correlation_mat = timepoint_correlation_mat
)
p_reference_overview <- p_reference_pca + p_reference_timepoint_correlation +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.08))

# 5.0 fit PCA and interpret loadings -----------------

timing_fit <- score_differentiation_timing(
  expression_matrix = vst,
  metadata = metadata,
  temporal_genes = temporal_genes,
  sample_id_col = 'sample_id',
  time_col = 'day_numeric'
)
pca_fit <- timing_fit$pca_fit
centroid_polyline <- timing_fit$centroid_polyline
pc1_day_correlation <- stats::cor(centroid_polyline$PC1, centroid_polyline$timepoint)
pc1_orientation <- if (is.finite(pc1_day_correlation) && pc1_day_correlation < 0) -1 else 1
pca_fit$x[, 'PC1'] <- pca_fit$x[, 'PC1'] * pc1_orientation
pca_fit$rotation[, 'PC1'] <- pca_fit$rotation[, 'PC1'] * pc1_orientation

timing_pca_variance_percent <- round(
  timing_fit$pca_variance$variance_percent[
    seq_len(timing_fit$n_pcs)
  ],
  1
)

timing_pca_data <- timing_fit$pca_coordinates
timing_pca_data$PC1 <- timing_pca_data$PC1 * pc1_orientation
timing_pca_data$day_numeric <- timing_pca_data$observed_time

early_cluster_labels <- utils::head(levels(gene_clusters), 2L)
late_cluster_labels <- utils::tail(levels(gene_clusters), 2L)
early_cluster_genes <- names(gene_clusters)[gene_clusters %in% early_cluster_labels]
late_cluster_genes <- names(gene_clusters)[gene_clusters %in% late_cluster_labels]

gene_sets <- list(
  All = temporal_genes,
  `C1+C2 (Early)` = early_cluster_genes,
  `C3+C4 (Late)` = late_cluster_genes
)

pca_results <- Map(calculate_pca, gene_sets, names(gene_sets))

p_all_pca <- make_pca_plot(pca_results[['All']], 'All')
p_early_pca <- make_pca_plot(pca_results[['C1+C2 (Early)']], 'C1+C2 (Early)')
p_late_pca <- make_pca_plot(pca_results[['C3+C4 (Late)']], 'C3+C4 (Late)')

p_all_pc1 <- make_pc1_time_plot(pca_results[['All']])
p_early_pc1 <- make_pc1_time_plot(pca_results[['C1+C2 (Early)']])
p_late_pc1 <- make_pc1_time_plot(pca_results[['C3+C4 (Late)']])

p_pca_day <- (
  p_all_pca | p_all_pc1
) / (
  p_early_pca | p_early_pc1
) / (
  p_late_pca | p_late_pc1
) +
  patchwork::plot_layout(widths = c(1.05, 0.95), heights = c(1, 1, 1)) &
  ggplot2::theme(plot.margin = ggplot2::margin(7, 8, 16, 8))

names(centroid_polyline)[names(centroid_polyline) == 'timepoint'] <- 'day_numeric'
centroid_polyline$PC1 <- centroid_polyline$PC1 * pc1_orientation

pc1_loadings <- pca_fit$rotation[, 'PC1']
top_pc1_pos_go <- names(
  sort(pc1_loadings, decreasing = TRUE)
)[seq_len(n_pc1_go_genes_per_direction)]
top_pc1_neg_go <- names(
  sort(pc1_loadings, decreasing = FALSE)
)[seq_len(n_pc1_go_genes_per_direction)]

pc1_loading_df <- make_loading_vector_data('PC1')
pc2_loading_df <- make_loading_vector_data('PC2')
p_pc1_loading_vectors <- make_loading_vector_plot(pc1_loading_df, 'PC1', 'Top PC1 loading genes')
p_pc2_loading_vectors <- make_loading_vector_plot(pc2_loading_df, 'PC2', 'Top PC2 loading genes')

pc1_go_pos <- run_go_enrichment(
  gene_ids = top_pc1_pos_go,
  universe_gene_ids = temporal_genes
)
pc1_go_neg <- run_go_enrichment(
  gene_ids = top_pc1_neg_go,
  universe_gene_ids = temporal_genes
)

p_pc1_go_pos <- make_go_dotplot(
  go_results = pc1_go_pos,
  plot_title = 'PC1+ GO:BP',
  point_color = pc1_positive_color
)
p_pc1_go_neg <- make_go_dotplot(
  go_results = pc1_go_neg,
  plot_title = 'PC1- GO:BP',
  point_color = pc1_negative_color
)

p_pc1_validation <- (p_pc1_loading_vectors | p_pc2_loading_vectors) / (p_pc1_go_neg | p_pc1_go_pos) +
  patchwork::plot_layout(heights = c(1.18, 0.68), widths = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    theme = ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      plot.tag = ggplot2::element_text(face = 'bold', family = figure_family, size = 8)
    )
  )

# 6.0 score samples on the centroid polyline -----------------

start_day <- min(days)
end_day <- max(days)
polyline_label_x_span <- diff(range(timing_pca_data$PC1, centroid_polyline$PC1))
polyline_label_y_span <- diff(range(timing_pca_data$PC2, centroid_polyline$PC2))

polyline_label_days <- c(1, 4, 7, 10, 15)
polyline_label_df <- centroid_polyline[
  centroid_polyline$day_numeric %in% polyline_label_days,
  c('day_numeric', 'PC1', 'PC2'),
  drop = FALSE
]
polyline_label_df$nudge_x <- c(-0.035, -0.055, 0.030, 0.035, 0.035) * polyline_label_x_span
polyline_label_df$nudge_y <- c(0.055, 0.045, 0.065, 0.040, -0.075) * polyline_label_y_span
polyline_label_df$label_x <- polyline_label_df$PC1 + polyline_label_df$nudge_x
polyline_label_df$label_y <- polyline_label_df$PC2 + polyline_label_df$nudge_y
polyline_label_df$day_label <- paste0('Day ', polyline_label_df$day_numeric)
polyline_endpoint_labels <- polyline_label_df[
  polyline_label_df$day_numeric %in% c(start_day, end_day),
  ,
  drop = FALSE
]
polyline_mid_labels <- polyline_label_df[
  !polyline_label_df$day_numeric %in% c(start_day, end_day),
  ,
  drop = FALSE
]

centroid_polyline_segments <- data.frame(
  PC1 = centroid_polyline$PC1[-nrow(centroid_polyline)],
  PC2 = centroid_polyline$PC2[-nrow(centroid_polyline)],
  PC1_end = centroid_polyline$PC1[-1L],
  PC2_end = centroid_polyline$PC2[-1L],
  stringsAsFactors = FALSE
)

timing_pca_data$differentiation_score <- timing_fit$scores$differentiation_score[
  match(timing_pca_data$sample_id, timing_fit$scores$sample_id)
]

p_timing_polyline <- ggplot2::ggplot(timing_pca_data, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = 0.95, alpha = 1) +
  ggplot2::geom_segment(
    data = centroid_polyline_segments,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      xend = PC1_end,
      yend = PC2_end
    ),
    inherit.aes = FALSE,
    color = 'black',
    linewidth = 0.48,
    arrow = grid::arrow(length = grid::unit(0.052, 'in'), type = 'closed')
  ) +
  ggplot2::geom_point(
    data = centroid_polyline,
    ggplot2::aes(PC1, PC2),
    inherit.aes = FALSE,
    color = 'black',
    size = 0.78
  ) +
  ggplot2::geom_label(
    data = polyline_mid_labels,
    ggplot2::aes(label_x, label_y, label = day_label),
    inherit.aes = FALSE,
    family = figure_family,
    fontface = 'bold',
    size = text_size(6.4),
    fill = 'white',
    color = 'black',
    linewidth = 0.14,
    label.padding = grid::unit(0.14, 'lines'),
    label.r = grid::unit(0.05, 'lines')
  ) +
  ggplot2::geom_label(
    data = polyline_endpoint_labels,
    ggplot2::aes(label_x, label_y, label = day_label),
    inherit.aes = FALSE,
    family = figure_family,
    fontface = 'bold',
    size = text_size(7.6),
    fill = 'white',
    color = 'black',
    linewidth = 0.14,
    label.padding = grid::unit(0.16, 'lines'),
    label.r = grid::unit(0.05, 'lines')
  ) +
  scale_color_by_day(name = 'Day') +
  ggplot2::labs(
    x = paste0('PC1 (', timing_pca_variance_percent[1], '%)'),
    y = paste0('PC2 (', timing_pca_variance_percent[2], '%)')
  ) +
  theme_tutorial() +
  theme_legend_overlay()

score_summary <- stats::aggregate(
  differentiation_score ~ day_numeric,
  data = timing_pca_data,
  FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
)

score_summary <- data.frame(
  day_numeric = score_summary$day_numeric,
  mean_score = score_summary$differentiation_score[, 'mean'],
  sd_score = score_summary$differentiation_score[, 'sd'],
  row.names = NULL
)
score_summary$lower_score <- score_summary$mean_score - score_summary$sd_score
score_summary$upper_score <- score_summary$mean_score + score_summary$sd_score
score_curve <- data.frame(
  day_numeric = seq(min(days), max(days), length.out = 220),
  mean_score = smooth_score_series(score_summary$mean_score),
  lower_score = smooth_score_series(score_summary$lower_score),
  upper_score = smooth_score_series(score_summary$upper_score)
)
score_curve_lower <- pmin(score_curve$lower_score, score_curve$upper_score)
score_curve_upper <- pmax(score_curve$lower_score, score_curve$upper_score)
score_curve$lower_score <- score_curve_lower
score_curve$upper_score <- score_curve_upper

p_score_by_day <- ggplot2::ggplot(timing_pca_data, ggplot2::aes(day_numeric, differentiation_score)) +
  ggplot2::geom_hline(
    yintercept = c(0, 1),
    linetype = 'dashed',
    linewidth = 0.24,
    color = 'black'
  ) +
  ggplot2::geom_ribbon(
    data = score_curve,
    ggplot2::aes(
      x = day_numeric,
      ymin = lower_score,
      ymax = upper_score
    ),
    inherit.aes = FALSE,
    fill = 'grey55',
    alpha = 0.24
  ) +
  ggplot2::geom_line(
    data = score_curve,
    ggplot2::aes(day_numeric, mean_score),
    inherit.aes = FALSE,
    linewidth = 0.52,
    color = 'black'
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = day_numeric),
    size = 0.95,
    alpha = 1
  ) +
  ggplot2::scale_x_continuous(breaks = days) +
  scale_color_by_day(name = 'Day') +
  ggplot2::labs(x = 'Differentiation day', y = 'Differentiation timing score') +
  theme_tutorial() +
  ggplot2::theme(
    legend.position = c(0.36, 0.78),
    legend.justification = c(0.5, 0.5),
    legend.direction = 'horizontal',
    legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    legend.key.width = grid::unit(0.30, 'in'),
    legend.title = ggplot2::element_text(
      family = figure_family,
      face = 'bold',
      size = 6,
      hjust = 0.5
    )
  )

# 7.0 load or run cell-line cross-validation -----------------

loo_validation <- NULL
if (file.exists(loo_validation_cache_path)) {
  expected_validation_settings <- .leave_one_cell_line_out_validation_settings(
    expression_cpm_cutoff,
    lrt_padj_cutoff,
    vst_dynamic_range_cutoff
  )
  cache_candidate <- tryCatch(
    readRDS(loo_validation_cache_path),
    error = function(error) NULL
  )
  cache_has_expected_structure <-
    is.list(cache_candidate) &&
    is.data.frame(cache_candidate$scores) &&
    is.data.frame(cache_candidate$summary) &&
    all(
      c('sample_id', 'heldout_cell_line', 'gene_set', 'predicted_day') %in%
        names(cache_candidate$scores)
    ) &&
    all(c('n_pcs', 'retained_variance_percent') %in% names(cache_candidate$summary))
  cache_matches_settings <- isTRUE(cache_has_expected_structure) &&
    identical(cache_candidate$settings, expected_validation_settings)
  cache_matches_samples <- isTRUE(cache_has_expected_structure) &&
    setequal(
      unique(cache_candidate$scores$sample_id),
      metadata$sample_id
    ) &&
    setequal(
      unique(as.character(cache_candidate$scores$heldout_cell_line)),
      unique(as.character(metadata$cell_line))
    )
  cache_matches_pca_contract <- isTRUE(cache_has_expected_structure) &&
    isTRUE(all(cache_candidate$summary$n_pcs >= 1L)) &&
    isTRUE(all(cache_candidate$summary$retained_variance_percent >= 99))

  if (cache_matches_settings && cache_matches_samples && cache_matches_pca_contract) {
    loo_validation <- cache_candidate
    message('Using optional validation cache: ', loo_validation_cache_path)
  } else {
    message('Ignoring incompatible validation cache: ', loo_validation_cache_path)
  }
}
if (is.null(loo_validation)) {
  message('Validation cache unavailable; running leave-one-cell-line-out validation.')
  loo_validation <- run_leave_one_cell_line_out_validation(
    counts = counts,
    vst = vst,
    metadata = metadata,
    expression_cpm_cutoff = expression_cpm_cutoff,
    lrt_padj_cutoff = lrt_padj_cutoff,
    vst_dynamic_range_cutoff = vst_dynamic_range_cutoff
  )
}
loo_all_temporal_scores <- loo_validation$scores[
  loo_validation$scores$gene_set == 'All temporal',
  ,
  drop = FALSE
]
loo_all_temporal_scores$heldout_cell_line <- factor(
  loo_all_temporal_scores$heldout_cell_line,
  levels = sort(unique(loo_all_temporal_scores$heldout_cell_line))
)
loo_display_excluded_cell_line <- '19190'
loo_display_scores <- loo_all_temporal_scores[
  !loo_all_temporal_scores$heldout_cell_line %in% loo_display_excluded_cell_line,
  ,
  drop = FALSE
]
loo_display_scores$heldout_cell_line <- factor(
  loo_display_scores$heldout_cell_line,
  levels = sort(unique(loo_display_scores$heldout_cell_line))
)
loo_display_cell_lines <- levels(loo_display_scores$heldout_cell_line)
loo_cell_line_ncol <- 4L
loo_cell_line_nrow <- ceiling(length(loo_display_cell_lines) / loo_cell_line_ncol)
loo_cell_line_plots <- lapply(seq_along(loo_display_cell_lines), function(i) {
  cell_line_id <- loo_display_cell_lines[[i]]
  cell_line_scores <- loo_display_scores[
    loo_display_scores$heldout_cell_line == cell_line_id,
    ,
    drop = FALSE
  ]
  cell_line_day_scores <- stats::aggregate(
    predicted_day ~ actual_day,
    data = cell_line_scores,
    FUN = mean
  )
  cell_line_day_counts <- stats::aggregate(
    predicted_day ~ actual_day,
    data = cell_line_scores,
    FUN = length
  )
  cell_line_day_scores$n_samples <- cell_line_day_counts$predicted_day
  row_index <- ceiling(i / loo_cell_line_ncol)
  col_index <- ((i - 1L) %% loo_cell_line_ncol) + 1L
  show_y_title <- col_index == 1L
  show_x_title <- row_index == loo_cell_line_nrow
  cell_line_r <- stats::cor(
    cell_line_scores$actual_day,
    cell_line_scores$predicted_day,
    use = 'complete.obs'
  )
  cell_line_r_label <- sprintf(
    'r = %.2f\nr\u00b2 = %.2f',
    cell_line_r,
    cell_line_r^2
  )

  ggplot2::ggplot(cell_line_scores, ggplot2::aes(actual_day, predicted_day)) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = 'black',
      linetype = 'dashed',
      linewidth = 0.24
    ) +
    ggplot2::geom_point(
      color = 'black',
      size = 0.55,
      alpha = 1
    ) +
    ggplot2::geom_smooth(
      data = cell_line_day_scores,
      mapping = ggplot2::aes(actual_day, predicted_day, weight = n_samples),
      inherit.aes = FALSE,
      method = 'loess',
      formula = y ~ x,
      se = FALSE,
      span = 0.78,
      color = '#1260A4',
      linewidth = 0.52
    ) +
    ggplot2::annotate(
      'text',
      x = min(days),
      y = max(days),
      label = cell_line_r_label,
      hjust = 0,
      vjust = 1,
      family = figure_family,
      size = text_size(5.1),
      lineheight = 0.86
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    ggplot2::scale_y_continuous(breaks = days) +
    ggplot2::coord_cartesian(xlim = range(days), ylim = range(days)) +
    ggplot2::labs(
      title = paste('Cell line', cell_line_id),
      x = if (show_x_title) 'Actual differentiation day' else NULL,
      y = if (show_y_title) 'Predicted differentiation day' else NULL
    ) +
    theme_tutorial() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = 'plain',
        hjust = 0,
        size = 7,
        family = figure_family,
        color = 'black'
      ),
      axis.text.x = ggplot2::element_text(size = 6, color = 'black', family = figure_family),
      axis.text.y = ggplot2::element_text(size = 6, color = 'black', family = figure_family),
      axis.title.x = ggplot2::element_text(size = 7, color = 'black', family = figure_family),
      axis.title.y = ggplot2::element_text(size = 7, color = 'black', family = figure_family),
      plot.margin = ggplot2::margin(5, 6, 6, 5)
    )
})
p_loo_cell_line_predictions <- patchwork::wrap_plots(
  loo_cell_line_plots,
  ncol = loo_cell_line_ncol
)

loo_r <- stats::cor(
  loo_all_temporal_scores$actual_day,
  loo_all_temporal_scores$predicted_day,
  use = 'complete.obs'
)
loo_abs_error <- abs(loo_all_temporal_scores$residual)
loo_accuracy_summary <- data.frame(
  metric = c(
    'Held-out cell lines',
    'Held-out samples',
    'Correlation between actual and predicted day',
    'Squared Pearson correlation',
    'Mean absolute error',
    'Median absolute error',
    'Predictions within 1 day',
    'Predictions within 2 days'
  ),
  value = c(
    format(length(unique(loo_all_temporal_scores$heldout_cell_line)), big.mark = ','),
    format(nrow(loo_all_temporal_scores), big.mark = ','),
    sprintf('%.3f', loo_r),
    sprintf('%.3f', loo_r^2),
    sprintf('%.2f days', mean(loo_abs_error, na.rm = TRUE)),
    sprintf('%.2f days', stats::median(loo_abs_error, na.rm = TRUE)),
    sprintf('%.1f%%', mean(loo_abs_error <= 1, na.rm = TRUE) * 100),
    sprintf('%.1f%%', mean(loo_abs_error <= 2, na.rm = TRUE) * 100)
  ),
  stringsAsFactors = FALSE
)

loo_regression_label <- format_regression_label(
  x = loo_all_temporal_scores$actual_day,
  y = loo_all_temporal_scores$predicted_day
)

p_loo_summary <- ggplot2::ggplot(
  loo_all_temporal_scores,
  ggplot2::aes(actual_day, predicted_day)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    color = 'black',
    linetype = 'dashed',
    linewidth = 0.24
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = heldout_cell_line),
    size = 0.95,
    alpha = 1
  ) +
  ggplot2::geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = FALSE,
    color = 'black',
    linewidth = 0.45
  ) +
  ggplot2::annotate(
    'text',
    x = min(loo_all_temporal_scores$actual_day),
    y = max(loo_all_temporal_scores$predicted_day),
    label = loo_regression_label,
    hjust = 0,
    vjust = 1,
    family = figure_family,
    size = text_size(5.5),
    lineheight = 0.95
  ) +
  ggplot2::scale_color_viridis_d(name = 'Held-out cell line', option = 'D') +
  ggplot2::scale_x_continuous(breaks = days) +
  ggplot2::scale_y_continuous(breaks = days) +
  ggplot2::labs(x = 'Actual differentiation day', y = 'Predicted differentiation day') +
  theme_tutorial() +
  ggplot2::theme(legend.position = 'none')

loo_timepoint_errors <- loo_all_temporal_scores
loo_timepoint_errors$absolute_error <- abs(loo_timepoint_errors$residual)
loo_timepoint_errors$actual_day_factor <- factor(
  loo_timepoint_errors$actual_day,
  levels = days
)
loo_timepoint_summary <- base::do.call(
  rbind,
  lapply(days, function(day_value) {
    day_errors <- loo_timepoint_errors$absolute_error[
      loo_timepoint_errors$actual_day == day_value
    ]
    data.frame(
      actual_day_factor = factor(day_value, levels = days),
      mean_absolute_error = mean(day_errors),
      sd_absolute_error = stats::sd(day_errors)
    )
  })
)
loo_timepoint_summary$lower_error <- pmax(
  0,
  loo_timepoint_summary$mean_absolute_error - loo_timepoint_summary$sd_absolute_error
)
loo_timepoint_summary$upper_error <-
  loo_timepoint_summary$mean_absolute_error + loo_timepoint_summary$sd_absolute_error

p_loo_timepoint_accuracy <- ggplot2::ggplot(
  loo_timepoint_summary,
  ggplot2::aes(actual_day_factor, mean_absolute_error)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = 'black',
    linewidth = 0.22,
    linetype = 'dashed'
  ) +
  ggplot2::geom_col(
    width = 0.58,
    fill = '#5F89A9',
    color = '#294B63',
    linewidth = 0.28
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower_error, ymax = upper_error),
    width = 0.20,
    linewidth = 0.32,
    color = '#17212B'
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(
    x = 'Actual differentiation day',
    y = 'Absolute prediction error (days)'
  ) +
  ggplot2::theme_classic(base_size = 7, base_family = figure_family) +
  ggplot2::theme(
    text = ggplot2::element_text(family = figure_family, color = 'black'),
    axis.title = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(size = 6),
    axis.text.y = ggplot2::element_text(size = 6),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )
# 8.0 save tutorial figures -----------------

save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_reference_pca_and_day_correlation'),
  plot = p_reference_overview,
  width = reference_overview_figure_width,
  height = reference_overview_figure_height
)
save_drawn_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_temporal_heatmap'),
  draw_fn = function() {
    ComplexHeatmap::ht_opt(
      ROW_ANNO_PADDING = grid::unit(0, 'mm'),
      DENDROGRAM_PADDING = grid::unit(0, 'mm'),
      TITLE_PADDING = grid::unit(c(1, 1), 'pt')
    )
    on.exit(ComplexHeatmap::ht_opt(RESET = TRUE), add = TRUE)
    ComplexHeatmap::draw(
      p_lrt_heatmap,
      heatmap_legend_side = 'right',
      annotation_legend_side = 'bottom'
    )
  },
  width = temporal_heatmap_figure_width,
  height = temporal_heatmap_figure_height,
  png_scale = 2
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_temporal_clusters'),
  plot = p_cluster_trajectories,
  width = temporal_clusters_figure_width,
  height = temporal_clusters_figure_height
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_pca_day'),
  plot = p_pca_day,
  width = pca_day_figure_width,
  height = pca_day_figure_height
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_pc1_validation'),
  plot = p_pc1_validation,
  width = pc1_validation_figure_width,
  height = pc1_validation_figure_height
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_timing_polyline'),
  plot = p_timing_polyline,
  width = timing_polyline_figure_width,
  height = timing_polyline_figure_height
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_score_by_day'),
  plot = p_score_by_day,
  width = score_by_day_figure_width,
  height = score_by_day_figure_height
)

save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_loo_cell_line_predictions'),
  plot = p_loo_cell_line_predictions,
  width = loo_cell_line_predictions_figure_width,
  height = loo_cell_line_predictions_figure_height
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_loo_summary'),
  plot = p_loo_summary,
  width = loo_summary_figure_width,
  height = loo_summary_figure_height
)
save_figure(
  path_stub = file.path(docs_figure_dir, 'GSE122380_loo_timepoint_accuracy'),
  plot = p_loo_timepoint_accuracy,
  width = loo_timepoint_accuracy_figure_width,
  height = loo_timepoint_accuracy_figure_height
)
