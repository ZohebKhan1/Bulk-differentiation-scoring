# ----
# author:
# - Zoheb Khan
#
# script path:
# - scripts/01_build_tutorial_objects.R
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
# - tmp/GSE122380_leave_one_line_out_validation.rds
#
# outputs:
# - tmp/GSE122380_temporal_cluster_go_k4.rds
# - docs/assets/figures/GSE122380_reference_pca_and_day_correlation.{png,svg}
# - docs/assets/figures/GSE122380_temporal_heatmap.{png,svg}
# - docs/assets/figures/GSE122380_temporal_clusters.{png,svg}
# - docs/assets/figures/GSE122380_pca_day.{png,svg}
# - docs/assets/figures/GSE122380_pc1_validation.{png,svg}
# - docs/assets/figures/GSE122380_timing_polyline.{png,svg}
# - docs/assets/figures/GSE122380_score_by_day.{png,svg}
# - docs/assets/figures/GSE122380_loo_line_predictions.{png,svg}
# - docs/assets/figures/GSE122380_loo_summary.{png,svg}
# - docs/assets/figures/GSE122380_loo_timepoint_accuracy.{png,svg}
# ----

# 0.0 source packages and dependencies -----------------

required_packages <- base::c(
  'DESeq2',
  'ComplexHeatmap',
  'AnnotationDbi',
  'circlize',
  'clusterProfiler',
  'edgeR',
  'ggrepel',
  'ggplot2',
  'org.Hs.eg.db',
  'patchwork',
  'ragg',
  'scales',
  'svglite',
  'systemfonts',
  'viridis',
  'yaml'
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
n_timing_pcs <- analysis_params$timing_score$n_pcs

docs_figure_dir <- 'docs/assets/figures'
tutorial_font_dir <- 'tutorial/assets/fonts'
loo_validation_cache_path <- 'tmp/GSE122380_leave_one_line_out_validation.rds'
cluster_go_cache_path <- 'tmp/GSE122380_temporal_cluster_go_k4.rds'

figure_dpi <- 600
figure_scale <- 1
figure_font_scale <- 1
figure_geom_scale <- 1
figure_family <- 'GSE122380 Nimbus Sans'
panel_tag_family <- figure_family
svg_figure_family <- 'Nimbus Sans'
fs <- function(size) size * figure_font_scale
fd <- function(size) size * figure_scale
gs <- function(size) size * figure_geom_scale
gfs <- function(size) fs(size) / ggplot2::.pt

n_heatmap_genes <- 1500L
n_temporal_clusters <- 4L
n_pc1_loading_plot_genes_per_direction <- 10L
n_pc1_go_genes_per_direction <- 500L
pca_gene_fraction <- 0.10
pc1_negative_color <- '#1E40AF'
pc1_positive_color <- '#A80000'
go_top_n_terms <- 5L
go_term_padj_cutoff <- 0.05
go_min_gene_count_in_term <- 10L
go_min_genes_in_go_db <- 26L
go_max_genes_in_go_db <- 499L
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

# 1.2 read direct inputs -----------------

GSE122380_data <- load_GSE122380_data()
meta <- GSE122380_data$metadata
counts <- GSE122380_data$counts
vst <- GSE122380_data$vst

# 1.3 define local helpers -----------------

theme_pub <- function(base_size = fs(7)) {
  ggplot2::theme_classic(base_size = base_size, base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = 'black', family = figure_family),
      axis.title = ggplot2::element_text(color = 'black', family = figure_family, size = fs(7)),
      axis.text = ggplot2::element_text(color = 'black', family = figure_family, size = fs(6)),
      axis.line = ggplot2::element_line(color = 'black', linewidth = gs(0.24)),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = gs(0.24)),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'bold',
        size = fs(6)
      ),
      legend.text = ggplot2::element_text(color = 'black', family = figure_family, size = fs(5.5)),
      strip.text = ggplot2::element_text(color = 'black', family = figure_family, size = fs(6)),
      plot.title = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'plain',
        size = fs(7)
      ),
      plot.subtitle = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        face = 'plain',
        size = fs(6)
      ),
      plot.tag = ggplot2::element_text(
        color = 'black',
        family = panel_tag_family,
        face = 'bold',
        size = fs(8)
      )
    )
}

day_color_scale <- function(name = 'Day', option = 'D') {
  ggplot2::scale_color_viridis_c(
    name = name,
    option = option,
    breaks = base::c(1, 5, 10, 15),
    guide = ggplot2::guide_colorbar(
      direction = 'horizontal',
      title.position = 'top',
      title.hjust = 0.5,
      barwidth = grid::unit(0.78, 'in'),
      barheight = grid::unit(0.06, 'in'),
      frame.colour = 'black',
      frame.linewidth = gs(0.18),
      ticks = FALSE,
      theme = ggplot2::theme(
        legend.title = ggplot2::element_text(size = fs(6), face = 'bold', family = figure_family),
        legend.text = ggplot2::element_text(size = fs(5.5), family = figure_family),
        legend.ticks = ggplot2::element_blank(),
        legend.ticks.length = grid::unit(0, 'pt')
      )
    )
  )
}

legend_overlay_theme <- function() {
  ggplot2::theme(
    legend.position = base::c(0.98, 0.98),
    legend.justification = base::c(1, 1),
    legend.direction = 'horizontal',
    legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key.width = grid::unit(0.30, 'in')
  )
}

legend_top_theme <- function() {
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

save_figure <- function(path, plot, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    dpi = figure_dpi,
    bg = 'white',
    limitsize = FALSE
  )
}

svg_font_faces <- function() {
  list(
    svglite::font_face(
      svg_figure_family,
      otf = '../fonts/NimbusSans-Regular.otf',
      weight = 400
    ),
    svglite::font_face(
      svg_figure_family,
      otf = '../fonts/NimbusSans-Bold.otf',
      weight = 700
    ),
    svglite::font_face(
      svg_figure_family,
      otf = '../fonts/NimbusSans-Italic.otf',
      style = 'italic',
      weight = 400
    ),
    svglite::font_face(
      svg_figure_family,
      otf = '../fonts/NimbusSans-BoldItalic.otf',
      style = 'italic',
      weight = 700
    )
  )
}

save_svg <- function(path, plot, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = svglite::svglite,
    web_fonts = svg_font_faces(),
    fix_text_size = FALSE,
    bg = 'white',
    limitsize = FALSE
  )
}

save_figure_pair <- function(path_stub, plot, width, height) {
  width <- fd(width)
  height <- fd(height)
  save_figure(
    path = paste0(path_stub, '.png'),
    plot = plot,
    width = width,
    height = height
  )
  save_svg(
    path = paste0(path_stub, '.svg'),
    plot = plot,
    width = width,
    height = height
  )
}

save_drawn_png <- function(path, draw_fn, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(
    filename = path,
    width = width,
    height = height,
    units = 'in',
    res = figure_dpi,
    background = 'white'
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  draw_fn()
}

save_drawn_svg <- function(path, draw_fn, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  svglite::svglite(
    file = path,
    width = width,
    height = height,
    bg = 'white',
    web_fonts = svg_font_faces(),
    fix_text_size = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  draw_fn()
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

save_drawn_figure_pair <- function(path_stub, draw_fn, width, height, png_scale = 1) {
  width <- fd(width)
  height <- fd(height)
  save_drawn_png(
    path = paste0(path_stub, '.png'),
    draw_fn = draw_fn,
    width = width * png_scale,
    height = height * png_scale
  )
  save_drawn_svg(
    path = paste0(path_stub, '.svg'),
    draw_fn = draw_fn,
    width = width,
    height = height
  )
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

orient_pca_fit_by_day <- function(pca_fit, sample_days) {
  pc1_day_cor <- stats::cor(
    pca_fit$x[, 'PC1'],
    sample_days[base::match(base::rownames(pca_fit$x), base::names(sample_days))],
    use = 'pairwise.complete.obs'
  )
  if (base::is.finite(pc1_day_cor) && pc1_day_cor < 0) {
    pca_fit$x[, 'PC1'] <- -pca_fit$x[, 'PC1']
    pca_fit$rotation[, 'PC1'] <- -pca_fit$rotation[, 'PC1']
  }
  pca_fit
}

get_pca_color_option <- function(set_label) {
  if (base::identical(set_label, 'C1+C2 (Early)')) {
    return('inferno')
  }
  if (base::identical(set_label, 'C3+C4 (Late)')) {
    return('mako')
  }
  'viridis'
}

format_lm_label <- function(x, y) {
  r_value <- stats::cor(x, y, use = 'pairwise.complete.obs')
  base::sprintf("r = %.2f\nr\u00b2 = %.2f", r_value, r_value^2)
}

format_go_label <- function(x, width = 34L) {
  base::vapply(
    x,
    function(term) base::paste(base::strwrap(term, width = width), collapse = '\n'),
    base::character(1)
  )
}

make_go_count_legend_breaks <- function(counts) {
  counts <- counts[base::is.finite(counts)]
  if (base::length(counts) == 0L) {
    return(base::numeric(0))
  }
  as.numeric(stats::quantile(
    counts,
    probs = base::c(0, 0.25, 0.75, 1),
    type = 7,
    names = FALSE
  ))
}

format_go_count_labels <- function(counts) {
  base::ifelse(
    base::abs(counts - base::round(counts)) < 0.01,
    base::as.character(base::round(counts)),
    base::formatC(counts, format = 'f', digits = 1)
  )
}

run_go_enrichment <- function(gene_ids, universe_gene_ids, set_label) {
  clean_gene_ids <- base::unique(base::sub('\\..*$', '', gene_ids))
  clean_universe <- base::unique(base::sub('\\..*$', '', universe_gene_ids))

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

  go_results <- base::as.data.frame(go_fit)
  if (base::nrow(go_results) == 0L) {
    return(go_results)
  }

  go_results <- go_results[
    !base::is.na(go_results$p.adjust) &
      go_results$p.adjust < go_term_padj_cutoff &
      go_results$Count >= go_min_gene_count_in_term,
    ,
    drop = FALSE
  ]
  if (base::nrow(go_results) == 0L) {
    return(go_results)
  }

  go_results$set <- set_label
  go_results$neg_log10_pvalue <- -base::log10(go_results$pvalue)
  go_results <- go_results[
    base::order(go_results$neg_log10_pvalue, go_results$Count, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  utils::head(go_results, go_top_n_terms)
}

make_go_dotplot <- function(go_results, plot_title, point_color, show_count_legend = TRUE) {
  if (base::nrow(go_results) == 0L) {
    empty_df <- base::data.frame(x = 0, y = plot_title)
    return(
      ggplot2::ggplot(empty_df, ggplot2::aes(x, y)) +
        ggplot2::geom_blank() +
        ggplot2::annotate(
          'text',
          x = 0,
          y = plot_title,
          label = 'No enriched GO terms',
          family = figure_family,
          size = gfs(5.5)
        ) +
        ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
        theme_pub() +
        ggplot2::theme(
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(
            color = point_color,
            face = 'bold',
            family = figure_family,
            size = fs(7)
          )
        )
    )
  }

  plot_df <- go_results
  plot_df$description_label <- format_go_label(plot_df$Description)
  plot_df$description_label <- base::factor(
    plot_df$description_label,
    levels = base::rev(plot_df$description_label)
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
      range = gs(base::c(1.1, 3.0)),
      guide = ggplot2::guide_legend(
        direction = 'vertical',
        title.position = 'top',
        title.hjust = 0,
        override.aes = base::list(alpha = 0.96)
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      x = '-log10(p)',
      y = NULL
    ) +
    theme_pub() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = point_color,
        face = 'bold',
        family = figure_family,
        size = fs(7)
      ),
      axis.text.y = ggplot2::element_text(
        color = 'black',
        family = figure_family,
        lineheight = 0.9,
        size = fs(5.5)
      ),
      legend.position = if (base::isTRUE(show_count_legend)) 'right' else 'none',
      legend.direction = 'vertical',
      legend.title = ggplot2::element_text(face = 'bold', family = figure_family, size = fs(6)),
      legend.text = ggplot2::element_text(family = figure_family, size = fs(5.5)),
      legend.key.size = grid::unit(0.26, 'cm'),
      legend.margin = ggplot2::margin(0, 0, 0, 2),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(4, 3, 4, 3)
    )
}

format_lm_equation_label <- function(x, y) {
  lm_fit <- stats::lm(y ~ x)
  coefs <- stats::coef(lm_fit)
  r_value <- stats::cor(x, y, use = 'pairwise.complete.obs')
  intercept <- coefs[[1]]
  slope <- coefs[[2]]
  sign_text <- if (slope < 0) '-' else '+'
  base::sprintf(
    "r = %.2f\nr\u00b2 = %.2f\ny = %.2f %s %.2fx",
    r_value,
    r_value^2,
    intercept,
    sign_text,
    base::abs(slope)
  )
}

plot_reference_pca <- function(pca_data, day_path, x_label, y_label, pca_gene_count) {
  reference_day_colors <- stats::setNames(reference_day_palette, base::seq_len(15))
  phase_colors <- base::c(
    'Pluripotent' = reference_day_colors[['1']],
    'Mesoderm' = reference_day_colors[['4']],
    'Immature cardiomyocyte' = reference_day_colors[['8']],
    'Mature cardiomyocyte' = reference_day_colors[['15']]
  )
  phase_label_fills <- phase_colors
  phase_label_fills[['Mature cardiomyocyte']] <- '#B8860B'

  pca_data$phase <- base::ifelse(
    pca_data$day_numeric <= 2,
    'Pluripotent',
    base::ifelse(
      pca_data$day_numeric <= 5,
      'Mesoderm',
      base::ifelse(
        pca_data$day_numeric <= 9,
        'Immature cardiomyocyte',
        'Mature cardiomyocyte'
      )
    )
  )
  pca_data$phase <- base::factor(pca_data$phase, levels = base::names(phase_colors))

  ellipse_df <- base::do.call(rbind, base::lapply(base::names(phase_colors), function(phase_name) {
    phase_data <- pca_data[pca_data$phase == phase_name, base::c('PC1', 'PC2'), drop = FALSE]
    if (base::nrow(phase_data) < 3L) {
      return(NULL)
    }

    center <- base::colMeans(phase_data)
    covariance <- stats::cov(phase_data)
    ellipse_angle <- base::seq(0, 2 * base::pi, length.out = 160)
    ellipse_circle <- base::rbind(base::cos(ellipse_angle), base::sin(ellipse_angle))
    ellipse_coords <- tryCatch(
      base::t(center + base::sqrt(stats::qchisq(0.80, df = 2)) *
        base::t(base::chol(covariance)) %*% ellipse_circle),
      error = function(e) NULL
    )
    if (base::is.null(ellipse_coords)) {
      return(NULL)
    }

    base::data.frame(
      phase = phase_name,
      PC1 = ellipse_coords[, 1],
      PC2 = ellipse_coords[, 2],
      stringsAsFactors = FALSE
    )
  }))

  ellipse_layers <- base::lapply(base::names(phase_colors), function(phase_name) {
    ggplot2::geom_path(
      data = ellipse_df[ellipse_df$phase == phase_name, , drop = FALSE],
      ggplot2::aes(PC1, PC2),
      inherit.aes = FALSE,
      color = phase_colors[[phase_name]],
      linewidth = gs(0.50),
      alpha = 0.9
    )
  })

  day_label_offsets <- base::data.frame(
    day_numeric = base::seq_len(15),
    nudge_x = base::c(-2, -4, -2, 0, 0, 0, -2, -2, -8, 5, -10, 8, -9, 10, 0),
    nudge_y = base::c(-6, 5, -5, 5, -5, 6, -7, 6, -5, 5, -9, 0, -2, -5, -12)
  )
  day_path <- base::merge(day_path, day_label_offsets, by = 'day_numeric', all.x = TRUE)
  day_path <- day_path[base::order(day_path$day_numeric), ]
  day_path$day_label <- base::paste0('D', day_path$day_numeric)
  day_path$label_x <- day_path$PC1 + day_path$nudge_x
  day_path$label_y <- day_path$PC2 + day_path$nudge_y

  phase_labels <- base::data.frame(
    day_numeric = base::c(1, 4, 8, 14),
    phase = base::c('Pluripotent', 'Mesoderm', 'Immature cardiomyocyte', 'Mature cardiomyocyte'),
    nudge_x = base::c(-4, -22, 8, 2),
    nudge_y = base::c(-12, 4, 19, -16),
    stringsAsFactors = FALSE
  )
  phase_labels <- base::merge(
    phase_labels,
    day_path[, base::c('day_numeric', 'PC1', 'PC2')],
    by = 'day_numeric',
    all.x = TRUE
  )
  phase_labels$label_x <- phase_labels$PC1 + phase_labels$nudge_x
  phase_labels$label_y <- phase_labels$PC2 + phase_labels$nudge_y
  phase_labels$phase <- base::factor(phase_labels$phase, levels = base::names(phase_colors))

  plot_bounds <- base::rbind(
    pca_data[, base::c('PC1', 'PC2')],
    ellipse_df[, base::c('PC1', 'PC2')],
    stats::setNames(day_path[, base::c('label_x', 'label_y')], base::c('PC1', 'PC2')),
    stats::setNames(phase_labels[, base::c('label_x', 'label_y')], base::c('PC1', 'PC2'))
  )
  pca_xlim <- base::range(plot_bounds$PC1, na.rm = TRUE) + base::c(-8, 8)
  pca_ylim <- base::range(plot_bounds$PC2, na.rm = TRUE) + base::c(-8, 8)

  ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2)) +
    ellipse_layers +
    ggplot2::geom_point(
      ggplot2::aes(color = day_numeric),
      size = gs(0.95),
      alpha = 0.92
    ) +
    ggplot2::geom_text(
      data = day_path,
      ggplot2::aes(label_x, label_y, label = day_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5),
      color = 'black'
    ) +
    ggplot2::geom_label(
      data = phase_labels,
      ggplot2::aes(label_x, label_y, label = phase, fill = phase),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5),
      color = 'white',
      linewidth = 0,
      label.r = grid::unit(0.08, 'lines'),
      label.padding = grid::unit(0.18, 'lines')
    ) +
    ggplot2::scale_color_gradientn(
      colors = reference_day_palette,
      name = 'Differentiation day',
      breaks = base::c(1, 5, 10, 15),
      guide = ggplot2::guide_colorbar(
        title.position = 'top',
        title.hjust = 0,
        barwidth = grid::unit(0.78, 'in'),
        barheight = grid::unit(0.06, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.18),
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
      subtitle = base::paste0('Top 10% variable genes (n=', scales::comma(pca_gene_count), ')'),
      tag = 'a',
      x = x_label,
      y = y_label
    ) +
    ggplot2::coord_cartesian(xlim = pca_xlim, ylim = pca_ylim, clip = 'off') +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      plot.title = ggplot2::element_text(family = figure_family, face = 'plain', size = fs(7)),
      plot.subtitle = ggplot2::element_text(family = figure_family, face = 'plain', size = fs(6)),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.title = ggplot2::element_text(family = figure_family, size = fs(7)),
      axis.text = ggplot2::element_text(family = figure_family, size = fs(6)),
      legend.position = base::c(0.98, 0.98),
      legend.justification = base::c(1, 1),
      legend.direction = 'horizontal',
      legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
      legend.margin = ggplot2::margin(1, 1, 1, 1),
      legend.title = ggplot2::element_text(family = figure_family, size = fs(6), face = 'bold'),
      legend.text = ggplot2::element_text(family = figure_family, size = fs(5.5)),
      legend.key.width = grid::unit(0.30, 'in'),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      plot.margin = ggplot2::margin(10, 8, 6, 10)
    )
}

make_timepoint_correlation_matrix <- function(expression_mat, metadata) {
  day_order <- base::sort(base::unique(metadata$day_numeric))
  day_mean_mat <- base::sapply(day_order, function(day_value) {
    sample_ids <- metadata$sample_id[metadata$day_numeric == day_value]
    base::rowMeans(expression_mat[, sample_ids, drop = FALSE], na.rm = TRUE)
  })
  base::colnames(day_mean_mat) <- base::paste0('D', day_order)

  stats::cor(day_mean_mat, method = 'pearson', use = 'pairwise.complete.obs')
}

get_day_colors <- function(days_to_color) {
  day_index <- base::round(scales::rescale(
    days_to_color,
    to = base::c(1, 256),
    from = base::range(days_to_color)
  ))
  stats::setNames(annotation_day_palette[day_index], base::paste0('D', days_to_color))
}

plot_timepoint_correlation_heatmap <- function(correlation_mat) {
  ordered_labels <- base::paste0('D', base::seq_len(base::ncol(correlation_mat)))
  y_labels <- base::rev(ordered_labels)
  n_labels <- base::length(ordered_labels)
  gap_after <- base::ceiling(n_labels / 2)
  gap_size <- 0.08
  add_heatmap_gap <- function(index) {
    index + base::ifelse(index > gap_after, gap_size, 0)
  }
  x_positions <- add_heatmap_gap(base::seq_along(ordered_labels))
  y_positions <- add_heatmap_gap(base::seq_along(y_labels))
  max_x_position <- base::max(x_positions)
  max_y_position <- base::max(y_positions)

  heatmap_data <- base::as.data.frame(base::as.table(correlation_mat[ordered_labels, y_labels]))
  base::names(heatmap_data) <- base::c('x_label', 'y_label', 'correlation')
  heatmap_data$x_index <- add_heatmap_gap(base::match(heatmap_data$x_label, ordered_labels))
  heatmap_data$y_index <- add_heatmap_gap(base::match(heatmap_data$y_label, y_labels))

  ordered_days <- base::as.integer(base::sub('^D', '', ordered_labels))
  y_days <- base::as.integer(base::sub('^D', '', y_labels))
  day_colors <- get_day_colors(base::sort(base::unique(ordered_days)))
  top_annotation <- base::data.frame(
    x_index = x_positions,
    y_index = max_y_position + 0.82,
    fill = base::unname(day_colors[ordered_labels])
  )
  left_annotation <- base::data.frame(
    x_index = 0.15,
    y_index = y_positions,
    fill = base::unname(day_colors[base::paste0('D', y_days)])
  )
  top_labels <- base::data.frame(
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
      linewidth = gs(0.20),
      arrow = grid::arrow(length = grid::unit(0.055, 'in'), type = 'closed')
    ) +
    ggplot2::annotate(
      'text',
      x = (max_x_position + 1) / 2,
      y = max_y_position + 2.84,
      label = 'Developmental time',
      family = figure_family,
      fontface = 'bold',
      size = gfs(5.5)
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
      size = gfs(5.5),
      angle = 45,
      hjust = 0,
      vjust = 0.5,
      color = 'grey25'
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      position = 'top',
      limits = base::c(-0.16, max_x_position + 0.51),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_positions,
      labels = y_labels,
      limits = base::c(0.49, max_y_position + 3.09),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_fill_gradientn(
      colors = correlation_palette,
      limits = base::c(0, 1),
      name = 'Pearson r',
      breaks = base::c(0, 0.5, 1),
      labels = function(x) {
        base::ifelse(
          x %% 1 == 0,
          base::formatC(x, format = 'f', digits = 0),
          base::formatC(x, format = 'f', digits = 1)
        )
      },
      guide = ggplot2::guide_colorbar(
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = grid::unit(1.25, 'in'),
        barheight = grid::unit(0.08, 'in'),
        frame.colour = 'black',
        frame.linewidth = gs(0.18),
        theme = ggplot2::theme(
          legend.ticks = ggplot2::element_blank(),
          legend.ticks.length = grid::unit(0, 'pt')
        )
      )
    ) +
    ggplot2::labs(title = NULL, tag = 'b', x = NULL, y = NULL) +
    ggplot2::coord_fixed(clip = 'off') +
    ggplot2::theme_minimal(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(family = panel_tag_family, face = 'bold', size = fs(8)),
      plot.tag.position = base::c(0.005, 0.995),
      axis.text.x = ggplot2::element_text(
        family = figure_family,
        size = fs(6),
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        margin = ggplot2::margin(b = -2)
      ),
      axis.text.y = ggplot2::element_text(family = figure_family, size = fs(6)),
      legend.position = 'bottom',
      legend.justification = 'center',
      legend.title = ggplot2::element_text(size = fs(6), face = 'bold', family = figure_family),
      legend.text = ggplot2::element_text(size = fs(5.5), family = figure_family),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-24, 0, 0, 34),
      legend.key.width = grid::unit(0.30, 'in'),
      legend.ticks = ggplot2::element_blank(),
      legend.ticks.length = grid::unit(0, 'pt'),
      plot.margin = ggplot2::margin(10, 8, 0, 10)
    )
}

# 2.0 select temporal genes -----------------

message('CODEX_STEP temporal_genes: selecting genes with CPM, LRT, and VST-range filters')

days <- sort(unique(meta$day_numeric))
sample_days <- stats::setNames(meta$day_numeric, meta$sample_id)
temporal_selection <- select_temporal_genes(
  counts = counts,
  vst = vst,
  metadata = meta,
  expression_cpm_cutoff = expression_cpm_cutoff,
  lrt_padj_cutoff = lrt_padj_cutoff,
  vst_dynamic_range_cutoff = vst_dynamic_range_cutoff
)
temporal_genes <- temporal_selection$temporal_genes
lrt_summary <- temporal_selection$summary

# 3.0 create temporal heatmap and clusters -----------------

message('CODEX_STEP heatmap: clustering temporal gene trajectories')

heatmap_genes <- utils::head(temporal_genes, base::min(n_heatmap_genes, base::length(temporal_genes)))
hm_mat <- vst[heatmap_genes, meta$sample_id, drop = FALSE]
hm_z <- base::t(base::scale(base::t(hm_mat)))
hm_z[!base::is.finite(hm_z)] <- 0
hm_z[hm_z > 2] <- 2
hm_z[hm_z < -2] <- -2

order_one_day <- function(day_value) {
  day_samples <- meta$sample_id[meta$day_numeric == day_value]
  if (base::length(day_samples) <= 2L) {
    return(day_samples)
  }
  day_dist <- stats::dist(base::t(hm_z[, day_samples, drop = FALSE]))
  day_samples[stats::hclust(day_dist, method = 'average')$order]
}

column_order <- base::unlist(base::lapply(days, order_one_day), use.names = FALSE)
hm_z <- hm_z[, column_order, drop = FALSE]
hm_meta <- meta[base::match(column_order, meta$sample_id), , drop = FALSE]

# k=4 is an interpretive sanity-check resolution. The score does not require
# clustering; clusters are used only to show that LRT-significant temporal genes
# span broad expression patterns and biological processes.
temporal_day_means <- base::t(base::vapply(heatmap_genes, function(gene) {
  base::vapply(days, function(day) {
    samples <- meta$sample_id[meta$day_numeric == day]
    base::mean(vst[gene, samples], na.rm = TRUE)
  }, base::numeric(1))
}, base::numeric(base::length(days))))
base::colnames(temporal_day_means) <- base::paste0('D', days)

temporal_day_z <- base::t(base::scale(base::t(temporal_day_means)))
temporal_day_z[!base::is.finite(temporal_day_z)] <- 0

smooth_temporal_trajectory <- function(gene_values) {
  smoothed <- tryCatch(
    {
      fit <- stats::loess(
        gene_values ~ days,
        span = 0.55,
        degree = 2,
        family = 'symmetric',
        control = stats::loess.control(surface = 'direct')
      )
      stats::predict(fit, newdata = days)
    },
    error = function(e) gene_values
  )
  smoothed[!base::is.finite(smoothed)] <- gene_values[!base::is.finite(smoothed)]
  smoothed
}

temporal_cluster_input <- base::t(base::apply(temporal_day_z, 1, smooth_temporal_trajectory))
temporal_cluster_input <- base::t(base::scale(base::t(temporal_cluster_input)))
temporal_cluster_input[!base::is.finite(temporal_cluster_input)] <- 0
temporal_cluster_tree <- stats::hclust(
  stats::dist(temporal_cluster_input),
  method = 'ward.D2'
)
raw_cluster <- stats::cutree(temporal_cluster_tree, k = n_temporal_clusters)

cluster_day_means <- base::do.call(rbind, base::lapply(base::sort(base::unique(raw_cluster)), function(cl) {
  genes <- base::names(raw_cluster)[raw_cluster == cl]
  day_values <- base::colMeans(temporal_day_z[genes, , drop = FALSE])
  base::data.frame(
    raw_cluster = cl,
    peak_day = days[base::which.max(day_values)],
    mean_signal = base::max(day_values),
    stringsAsFactors = FALSE
  )
}))

cluster_order <- cluster_day_means$raw_cluster[
  base::order(cluster_day_means$peak_day, -cluster_day_means$mean_signal)
]
cluster_map <- stats::setNames(base::paste0('C', base::seq_along(cluster_order)), cluster_order)
gene_cluster <- base::factor(
  cluster_map[base::as.character(raw_cluster)],
  levels = base::paste0('C', base::seq_along(cluster_order))
)
base::names(gene_cluster) <- base::names(raw_cluster)

cluster_colors <- stats::setNames(
  base::c('#0072B2', '#D55E00', '#009E73', '#CC79A7'),
  base::levels(gene_cluster)
)

col_fun <- circlize::colorRamp2(
  base::seq(-2, 2, length.out = 100),
  grDevices::colorRampPalette(correlation_palette)(100)
)

cluster_block_annotation <- ComplexHeatmap::rowAnnotation(
  cluster = ComplexHeatmap::anno_block(
    gp = grid::gpar(fill = base::unname(cluster_colors[base::levels(gene_cluster)]), col = NA),
    labels = base::paste('Cluster', base::seq_along(base::levels(gene_cluster))),
    labels_gp = grid::gpar(
      col = 'white',
      fontface = 'bold',
      fontfamily = figure_family,
      fontsize = fs(12)
    ),
    labels_rot = 270,
    width = grid::unit(12, 'mm'),
    show_name = FALSE
  ),
  width = grid::unit(12, 'mm'),
  show_annotation_name = FALSE
)

p_lrt_heatmap <- ComplexHeatmap::Heatmap(
  hm_z,
  name = 'Z-score',
  col = col_fun,
  right_annotation = cluster_block_annotation,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = gene_cluster[base::rownames(hm_z)],
  row_title = NULL,
  row_title_rot = 0,
  row_title_gp = grid::gpar(
    fontsize = fs(8),
    fontface = 'bold',
    fontfamily = figure_family,
    col = cluster_colors
  ),
  row_gap = grid::unit(1.8, 'mm'),
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = hm_meta$day_factor,
  column_gap = grid::unit(0, 'mm'),
  column_title = base::paste0('D', base::levels(hm_meta$day_factor)),
  column_title_gp = grid::gpar(
    fontsize = fs(15),
    fontface = 'bold',
    fontfamily = figure_family
  ),
  heatmap_legend_param = base::list(
    at = base::c(-2, 0, 2),
    labels = base::c('-2', '0', '2'),
    tick_length = grid::unit(0, 'mm'),
    border = 'black',
    legend_height = grid::unit(28, 'mm'),
    grid_width = grid::unit(4.8, 'mm'),
    legend_gp = grid::gpar(col = 'black', lwd = 0.5),
    title_gp = grid::gpar(fontface = 'bold', fontfamily = figure_family, fontsize = fs(10)),
    labels_gp = grid::gpar(fontfamily = figure_family, fontsize = fs(9))
  ),
  border = FALSE,
  rect_gp = grid::gpar(col = NA),
  use_raster = TRUE,
  raster_quality = 24
)

cluster_genes <- base::names(gene_cluster)
cluster_day_gene <- base::do.call(rbind, base::lapply(cluster_genes, function(gene) {
  base::data.frame(
    gene_id = gene,
    cluster = gene_cluster[[gene]],
    day_numeric = days,
    zscore = base::as.numeric(temporal_day_z[gene, ]),
    stringsAsFactors = FALSE
  )
}))

cluster_mean <- stats::aggregate(
  zscore ~ cluster + day_numeric,
  data = cluster_day_gene,
  FUN = base::mean
)

make_cluster_trajectory_plot <- function(cluster_name) {
  cluster_df <- cluster_day_gene[cluster_day_gene$cluster == cluster_name, ]
  mean_df <- cluster_mean[cluster_mean$cluster == cluster_name, ]
  cluster_title <- base::c(
    C1 = 'Cluster 1 (Stem cell/pluripotent)',
    C2 = 'Cluster 2 (Mesoderm/pluripotent)',
    C3 = 'Cluster 3 (Late cardiomyocyte maturation)',
    C4 = 'Cluster 4 (Early cardiomyocyte maturation)'
  )[[cluster_name]]
  y_limits <- base::list(
    C1 = base::c(-2, 3),
    C2 = base::c(-2, 2.5),
    C3 = base::c(-2, 2),
    C4 = base::c(-3, 2)
  )[[cluster_name]]
  y_breaks <- base::list(
    C1 = base::c(-2, 0, 3),
    C2 = base::c(-2, 2),
    C3 = base::c(-2, 0, 2),
    C4 = base::c(-3, 0, 2)
  )[[cluster_name]]

  ggplot2::ggplot(cluster_df, ggplot2::aes(day_numeric, zscore, group = gene_id)) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = 'black',
      linetype = 'dashed',
      linewidth = gs(0.24)
    ) +
    ggplot2::geom_line(color = 'gray72', linewidth = gs(0.18)) +
    ggplot2::geom_smooth(
      data = mean_df,
      ggplot2::aes(day_numeric, zscore),
      inherit.aes = FALSE,
      method = 'loess',
      formula = y ~ x,
      se = FALSE,
      span = 0.75,
      color = cluster_colors[[cluster_name]],
      linewidth = gs(0.75)
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(title = cluster_title, x = 'Differentiation day', y = 'Mean VST z-score') +
    theme_pub() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = cluster_colors[[cluster_name]],
        face = 'bold',
        family = figure_family,
        size = fs(7)
      )
    )
}

run_cluster_go_enrichment <- function(cluster_assignments, universe_gene_ids, cache_path) {
  cache_key <- list(
    input_md5 = GSE122380_data$input_md5,
    implementation_md5 = unname(tools::md5sum(c(
      'scripts/01_build_tutorial_objects.R',
      'config/analysis.yml'
    ))),
    package_versions = c(
      clusterProfiler = as.character(utils::packageVersion('clusterProfiler')),
      org.Hs.eg.db = as.character(utils::packageVersion('org.Hs.eg.db'))
    ),
    cluster_assignments = cluster_assignments,
    universe_gene_ids = sort(unique(sub('\\..*$', '', universe_gene_ids))),
    go_top_n_terms = go_top_n_terms,
    go_term_padj_cutoff = go_term_padj_cutoff,
    go_min_gene_count_in_term = go_min_gene_count_in_term,
    go_min_genes_in_go_db = go_min_genes_in_go_db,
    go_max_genes_in_go_db = go_max_genes_in_go_db
  )

  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (isTRUE(identical(cached$cache_key, cache_key))) {
      return(cached$cluster_go)
    }
  }

  message('CODEX_STEP go: running GO enrichment for temporal clusters')
  cluster_go <- lapply(levels(cluster_assignments), function(cluster_name) {
    run_go_enrichment(
      gene_ids = names(cluster_assignments)[cluster_assignments == cluster_name],
      universe_gene_ids = universe_gene_ids,
      set_label = cluster_name
    )
  })
  names(cluster_go) <- levels(cluster_assignments)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      cache_key = cache_key,
      cluster_go = cluster_go
    ),
    cache_path
  )
  cluster_go
}

cluster_go_results <- run_cluster_go_enrichment(
  cluster_assignments = gene_cluster,
  universe_gene_ids = heatmap_genes,
  cache_path = cluster_go_cache_path
)

make_cluster_go_plot <- function(cluster_name) {
  make_go_dotplot(
    go_results = cluster_go_results[[cluster_name]],
    plot_title = base::paste0(base::sub('^C', 'Cluster ', cluster_name), ' GO:BP'),
    point_color = cluster_colors[[cluster_name]],
    show_count_legend = TRUE
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = cluster_colors[[cluster_name]],
        face = 'bold',
        family = figure_family,
        size = fs(7)
      )
    )
}

cluster_panel_plots <- base::lapply(base::levels(gene_cluster), function(cluster_name) {
  (make_cluster_trajectory_plot(cluster_name) | make_cluster_go_plot(cluster_name)) +
    patchwork::plot_layout(widths = base::c(1.9, 0.72))
})

p_cluster_trajectories <- patchwork::wrap_plots(
  cluster_panel_plots,
  ncol = 1
) +
  patchwork::plot_annotation(tag_levels = 'a') &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = 'bold', family = panel_tag_family, size = fs(8)))

# 5.0 create reference dataset overview -----------------

pca_gene_vars <- base::apply(vst, 1, stats::var, na.rm = TRUE)
pca_gene_count <- base::ceiling(base::length(pca_gene_vars) * pca_gene_fraction)
pca_genes <- base::names(base::sort(pca_gene_vars, decreasing = TRUE))[base::seq_len(pca_gene_count)]
pca_input_mat <- vst[pca_genes, meta$sample_id, drop = FALSE]

overview_pca <- stats::prcomp(base::t(pca_input_mat), center = TRUE, scale. = FALSE)
overview_pca_df <- base::data.frame(
  sample_id = base::rownames(overview_pca$x),
  PC1 = overview_pca$x[, 1],
  PC2 = overview_pca$x[, 2],
  day_numeric = meta$day_numeric[base::match(base::rownames(overview_pca$x), meta$sample_id)],
  cell_line = meta$cell_line[base::match(base::rownames(overview_pca$x), meta$sample_id)],
  stringsAsFactors = FALSE
)

pc1_flip <- base::ifelse(
  stats::median(
    overview_pca_df$PC1[overview_pca_df$day_numeric == base::max(overview_pca_df$day_numeric)],
    na.rm = TRUE
  ) <
    stats::median(
      overview_pca_df$PC1[overview_pca_df$day_numeric == base::min(overview_pca_df$day_numeric)],
      na.rm = TRUE
    ),
  -1,
  1
)
overview_pca_df$PC1 <- overview_pca_df$PC1 * pc1_flip

overview_day_path <- stats::aggregate(
  overview_pca_df[, base::c('PC1', 'PC2')],
  by = base::list(day_numeric = overview_pca_df$day_numeric),
  FUN = stats::median,
  na.rm = TRUE
)
overview_day_path <- overview_day_path[base::order(overview_day_path$day_numeric), ]
overview_var_explained <- (overview_pca$sdev^2) / base::sum(overview_pca$sdev^2)
overview_x_label <- base::sprintf('PC1 (%.1f%%)', overview_var_explained[[1]] * 100)
overview_y_label <- base::sprintf('PC2 (%.1f%%)', overview_var_explained[[2]] * 100)

p_reference_pca <- plot_reference_pca(
  pca_data = overview_pca_df,
  day_path = overview_day_path,
  x_label = overview_x_label,
  y_label = overview_y_label,
  pca_gene_count = pca_gene_count
)

timepoint_correlation_mat <- make_timepoint_correlation_matrix(
  expression_mat = pca_input_mat,
  metadata = meta
)
p_reference_timepoint_correlation <- plot_timepoint_correlation_heatmap(
  correlation_mat = timepoint_correlation_mat
)
p_reference_overview <- p_reference_pca + p_reference_timepoint_correlation +
  patchwork::plot_layout(ncol = 2, widths = base::c(1, 1.08))

# 5.0 fit PCA and validate loadings -----------------

message('CODEX_STEP pca: fitting reference PCA and timing polyline')

timing_fit <- score_differentiation_timing(
  expression_matrix = vst,
  metadata = meta,
  temporal_genes = temporal_genes,
  sample_id_col = 'sample_id',
  time_col = 'day_numeric',
  n_pcs = n_timing_pcs
)
pca_fit_unoriented <- timing_fit$pca_fit
pca_fit <- pca_fit_unoriented
pca_fit <- orient_pca_fit_by_day(pca_fit, sample_days)
pc1_orientation <- if (isTRUE(all.equal(
  pca_fit$x[, 'PC1'],
  pca_fit_unoriented$x[, 'PC1']
))) {
  1
} else {
  -1
}

var_pct <- round(summary(pca_fit)$importance[2, seq_len(n_timing_pcs)] * 100, 1)

pca_df <- data.frame(
  sample_id = rownames(pca_fit$x),
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2],
  PC3 = pca_fit$x[, 3],
  meta[match(rownames(pca_fit$x), meta$sample_id), ],
  row.names = NULL
)

early_cluster_labels <- utils::head(base::levels(gene_cluster), 2L)
late_cluster_labels <- utils::tail(base::levels(gene_cluster), 2L)
early_cluster_genes <- base::names(gene_cluster)[gene_cluster %in% early_cluster_labels]
late_cluster_genes <- base::names(gene_cluster)[gene_cluster %in% late_cluster_labels]

gene_sets <- base::list(
  All = temporal_genes,
  `C1+C2 (Early)` = early_cluster_genes,
  `C3+C4 (Late)` = late_cluster_genes
)

base::stopifnot(base::all(base::vapply(gene_sets, base::length, base::integer(1)) > 1L))
base::stopifnot(base::all(base::unlist(gene_sets, use.names = FALSE) %in% base::rownames(vst)))

calculate_pca <- function(gene_ids, set_label) {
  set_pca_input <- base::t(vst[gene_ids, meta$sample_id, drop = FALSE])
  set_pca_fit <- stats::prcomp(set_pca_input, center = TRUE, scale. = FALSE)
  set_pca_fit <- orient_pca_fit_by_day(set_pca_fit, sample_days)
  set_pca_var <- base::round(base::summary(set_pca_fit)$importance[2, 1:3] * 100, 2)

  set_pca_df <- base::data.frame(
    sample_id = base::rownames(set_pca_fit$x),
    PC1 = set_pca_fit$x[, 1],
    PC2 = set_pca_fit$x[, 2],
    PC3 = set_pca_fit$x[, 3],
    day_numeric = meta$day_numeric[
      base::match(base::rownames(set_pca_fit$x), meta$sample_id)
    ],
    cell_line = meta$cell_line[
      base::match(base::rownames(set_pca_fit$x), meta$sample_id)
    ],
    gene_set = set_label,
    stringsAsFactors = FALSE
  )

  base::list(
    data = set_pca_df,
    variance_percent = set_pca_var,
    gene_count = base::length(gene_ids)
  )
}

pca_results <- base::Map(calculate_pca, gene_sets, base::names(gene_sets))

make_pca_plot <- function(pca_result, plot_title) {
  plot_pca_df <- pca_result$data
  plot_pca_var <- pca_result$variance_percent
  color_option <- get_pca_color_option(plot_title)

  ggplot2::ggplot(plot_pca_df, ggplot2::aes(PC1, PC2, color = day_numeric)) +
    ggplot2::geom_point(size = gs(0.95), alpha = 1) +
    day_color_scale(name = 'Day', option = color_option) +
    ggplot2::labs(
      title = plot_title,
      x = base::paste0('PC1 (', plot_pca_var[[1]], '%)'),
      y = base::paste0('PC2 (', plot_pca_var[[2]], '%)')
    ) +
    theme_pub() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = figure_family,
        face = 'plain',
        size = fs(7),
        margin = ggplot2::margin(b = -1)
      ),
      plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5)
    ) +
    legend_overlay_theme()
}

make_pc1_time_plot <- function(pca_result) {
  plot_pca_df <- pca_result$data
  color_option <- get_pca_color_option(plot_pca_df$gene_set[[1]])
  label_text <- format_lm_label(plot_pca_df$day_numeric, plot_pca_df$PC1)
  label_x <- base::min(plot_pca_df$day_numeric)
  label_y <- base::max(plot_pca_df$PC1, na.rm = TRUE)

  ggplot2::ggplot(plot_pca_df, ggplot2::aes(day_numeric, PC1)) +
    ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = gs(0.95), alpha = 1) +
    ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = gs(0.45)
    ) +
    ggplot2::annotate(
      'text',
      x = label_x,
      y = label_y,
      label = label_text,
      hjust = 0,
      vjust = 1,
      family = figure_family,
      size = gfs(7.5)
    ) +
    ggplot2::scale_x_continuous(breaks = days) +
    day_color_scale(name = 'Day', option = color_option) +
    ggplot2::labs(
      x = 'Differentiation day',
      y = 'PC1 position'
    ) +
    theme_pub() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 14)
    ) +
    legend_top_theme()
}

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
  patchwork::plot_layout(widths = base::c(1.05, 0.95), heights = base::c(1, 1, 1)) &
  ggplot2::theme(plot.margin = ggplot2::margin(7, 8, 16, 8))

centroid_polyline <- timing_fit$centroid_polyline
names(centroid_polyline)[names(centroid_polyline) == 'timepoint'] <- 'day_numeric'
centroid_polyline$PC1 <- centroid_polyline$PC1 * pc1_orientation

make_loading_vector_data <- function(component_name) {
  component_loadings <- pca_fit$rotation[, component_name]
  top_pos <- base::names(
    base::sort(component_loadings, decreasing = TRUE)
  )[base::seq_len(n_pc1_loading_plot_genes_per_direction)]
  top_neg <- base::names(
    base::sort(component_loadings, decreasing = FALSE)
  )[base::seq_len(n_pc1_loading_plot_genes_per_direction)]
  top_genes <- base::c(top_neg, top_pos)

  base::data.frame(
    gene_id = top_genes,
    gene_label = label_ensembl_genes(top_genes),
    PC1_loading = pca_fit$rotation[top_genes, 'PC1'],
    PC2_loading = pca_fit$rotation[top_genes, 'PC2'],
    direction = base::rep(
      base::c(base::paste0(component_name, '-'), base::paste0(component_name, '+')),
      each = n_pc1_loading_plot_genes_per_direction
    ),
    stringsAsFactors = FALSE
  )
}

make_loading_limit <- function(values, multiplier) {
  limit <- base::max(base::abs(values), na.rm = TRUE) * multiplier
  if (!base::is.finite(limit) || limit <= 0) {
    limit <- 0.01
  }
  base::c(-limit, limit)
}

make_loading_vector_plot <- function(loading_df, component_name, plot_title) {
  negative_label <- base::paste0(component_name, '-')
  positive_label <- base::paste0(component_name, '+')
  negative_df <- loading_df[loading_df$direction == negative_label, , drop = FALSE]
  positive_df <- loading_df[loading_df$direction == positive_label, , drop = FALSE]
  x_limits <- make_loading_limit(loading_df$PC1_loading, 1.45)
  y_limits <- make_loading_limit(loading_df$PC2_loading, 1.90)

  ggplot2::ggplot(loading_df) +
    ggplot2::geom_hline(yintercept = 0, color = 'grey65', linewidth = gs(0.20)) +
    ggplot2::geom_vline(xintercept = 0, color = 'grey65', linewidth = gs(0.20)) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = PC1_loading,
        yend = PC2_loading,
        color = direction
      ),
      linewidth = gs(0.24),
      arrow = grid::arrow(length = grid::unit(0.045, 'in'), type = 'closed')
    ) +
    ggrepel::geom_text_repel(
      data = negative_df,
      ggplot2::aes(PC1_loading, PC2_loading, label = gene_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'italic',
      size = gfs(5.7),
      color = pc1_negative_color,
      box.padding = grid::unit(0.08, 'lines'),
      point.padding = grid::unit(0.02, 'lines'),
      min.segment.length = 0,
      segment.color = NA,
      segment.linewidth = 0,
      max.overlaps = Inf,
      seed = 2026
    ) +
    ggrepel::geom_text_repel(
      data = positive_df,
      ggplot2::aes(PC1_loading, PC2_loading, label = gene_label),
      inherit.aes = FALSE,
      family = figure_family,
      fontface = 'italic',
      size = gfs(5.7),
      color = pc1_positive_color,
      box.padding = grid::unit(0.08, 'lines'),
      point.padding = grid::unit(0.02, 'lines'),
      min.segment.length = 0,
      segment.color = NA,
      segment.linewidth = 0,
      max.overlaps = Inf,
      seed = 2027
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(
        base::c(pc1_negative_color, pc1_positive_color),
        base::c(negative_label, positive_label)
      ),
      guide = 'none'
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'on') +
    ggplot2::labs(
      title = plot_title,
      x = 'PC1 loading',
      y = 'PC2 loading'
    ) +
    theme_pub() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = figure_family,
        face = 'bold',
        size = fs(7),
        hjust = 0.5
      ),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0)),
      plot.margin = ggplot2::margin(4, 3, 4, 3)
    )
}

pc1_loadings <- pca_fit$rotation[, 'PC1']
top_pc1_pos_go <- base::names(
  base::sort(pc1_loadings, decreasing = TRUE)
)[base::seq_len(n_pc1_go_genes_per_direction)]
top_pc1_neg_go <- base::names(
  base::sort(pc1_loadings, decreasing = FALSE)
)[base::seq_len(n_pc1_go_genes_per_direction)]

pc1_loading_df <- make_loading_vector_data('PC1')
pc2_loading_df <- make_loading_vector_data('PC2')
p_pc1_loading_vectors <- make_loading_vector_plot(pc1_loading_df, 'PC1', 'Top PC1 loading genes')
p_pc2_loading_vectors <- make_loading_vector_plot(pc2_loading_df, 'PC2', 'Top PC2 loading genes')

base::message('CODEX_STEP go: running GO enrichment for top PC1 loading genes')

pc1_go_pos <- run_go_enrichment(
  gene_ids = top_pc1_pos_go,
  universe_gene_ids = temporal_genes,
  set_label = 'PC1+'
)
pc1_go_neg <- run_go_enrichment(
  gene_ids = top_pc1_neg_go,
  universe_gene_ids = temporal_genes,
  set_label = 'PC1-'
)

p_pc1_go_pos <- make_go_dotplot(
  go_results = pc1_go_pos,
  plot_title = 'PC1+ GO:BP',
  point_color = pc1_positive_color,
  show_count_legend = TRUE
)
p_pc1_go_neg <- make_go_dotplot(
  go_results = pc1_go_neg,
  plot_title = 'PC1- GO:BP',
  point_color = pc1_negative_color,
  show_count_legend = TRUE
)

p_pc1_validation <- (p_pc1_loading_vectors | p_pc2_loading_vectors) / (p_pc1_go_neg | p_pc1_go_pos) +
  patchwork::plot_layout(heights = base::c(1.18, 0.68), widths = base::c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    theme = ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      plot.tag = ggplot2::element_text(face = 'bold', family = panel_tag_family, size = fs(8))
    )
  )

# 7.0 score samples on the centroid-polyline timing trajectory -----------------

start_day <- base::min(days)
end_day <- base::max(days)
polyline_label_x_span <- base::diff(base::range(pca_df$PC1, centroid_polyline$PC1))
polyline_label_y_span <- base::diff(base::range(pca_df$PC2, centroid_polyline$PC2))

polyline_label_days <- base::c(1, 4, 7, 10, 15)
polyline_label_df <- centroid_polyline[
  centroid_polyline$day_numeric %in% polyline_label_days,
  base::c('day_numeric', 'PC1', 'PC2'),
  drop = FALSE
]
polyline_label_df$nudge_x <- base::c(-0.035, -0.055, 0.030, 0.035, 0.035) * polyline_label_x_span
polyline_label_df$nudge_y <- base::c(0.055, 0.045, 0.065, 0.040, -0.075) * polyline_label_y_span
polyline_label_df$label_x <- polyline_label_df$PC1 + polyline_label_df$nudge_x
polyline_label_df$label_y <- polyline_label_df$PC2 + polyline_label_df$nudge_y
polyline_label_df$day_label <- base::paste0('Day ', polyline_label_df$day_numeric)
polyline_endpoint_labels <- polyline_label_df[
  polyline_label_df$day_numeric %in% base::c(start_day, end_day),
  ,
  drop = FALSE
]
polyline_mid_labels <- polyline_label_df[
  !polyline_label_df$day_numeric %in% base::c(start_day, end_day),
  ,
  drop = FALSE
]

centroid_polyline_segments <- base::data.frame(
  PC1 = centroid_polyline$PC1[-base::nrow(centroid_polyline)],
  PC2 = centroid_polyline$PC2[-base::nrow(centroid_polyline)],
  PC1_end = centroid_polyline$PC1[-1L],
  PC2_end = centroid_polyline$PC2[-1L],
  stringsAsFactors = FALSE
)

score_rows <- match(pca_df$sample_id, timing_fit$scores$sample_id)
stopifnot(!anyNA(score_rows))
pca_df$differentiation_score <- timing_fit$scores$differentiation_score[score_rows]
pca_df$predicted_day <- timing_fit$scores$predicted_time[score_rows]
pca_df$polyline_segment_start_day <- timing_fit$scores$nearest_segment_start[score_rows]
pca_df$polyline_segment_end_day <- timing_fit$scores$nearest_segment_end[score_rows]
pca_df$polyline_segment_fraction <- timing_fit$scores$segment_fraction[score_rows]
pca_df$polyline_squared_distance <- timing_fit$scores$squared_distance[score_rows]

p_timing_polyline <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2)) +
  ggplot2::geom_point(ggplot2::aes(color = day_numeric), size = gs(0.95), alpha = 1) +
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
    linewidth = gs(0.48),
    arrow = grid::arrow(length = grid::unit(0.052, 'in'), type = 'closed')
  ) +
  ggplot2::geom_point(
    data = centroid_polyline,
    ggplot2::aes(PC1, PC2),
    inherit.aes = FALSE,
    color = 'black',
    size = gs(0.78)
  ) +
  ggplot2::geom_label(
    data = polyline_mid_labels,
    ggplot2::aes(label_x, label_y, label = day_label),
    inherit.aes = FALSE,
    family = figure_family,
    fontface = 'bold',
    size = gfs(6.4),
    fill = 'white',
    color = 'black',
    linewidth = gs(0.14),
    label.padding = grid::unit(0.14, 'lines'),
    label.r = grid::unit(0.05, 'lines')
  ) +
  ggplot2::geom_label(
    data = polyline_endpoint_labels,
    ggplot2::aes(label_x, label_y, label = day_label),
    inherit.aes = FALSE,
    family = figure_family,
    fontface = 'bold',
    size = gfs(7.6),
    fill = 'white',
    color = 'black',
    linewidth = gs(0.14),
    label.padding = grid::unit(0.16, 'lines'),
    label.r = grid::unit(0.05, 'lines')
  ) +
  day_color_scale(name = 'Day') +
  ggplot2::labs(
    x = base::paste0("PC1 (", var_pct[1], "%)"),
    y = base::paste0("PC2 (", var_pct[2], "%)")
  ) +
  theme_pub() +
  legend_overlay_theme()

score_summary <- stats::aggregate(
  differentiation_score ~ day_numeric,
  data = pca_df,
  FUN = function(x) base::c(mean = base::mean(x), sd = stats::sd(x), n = base::length(x))
)

score_summary <- base::data.frame(
  day_numeric = score_summary$day_numeric,
  mean_score = score_summary$differentiation_score[, "mean"],
  sd_score = score_summary$differentiation_score[, "sd"],
  n = score_summary$differentiation_score[, "n"],
  row.names = NULL
)
score_summary$lower_score <- score_summary$mean_score - score_summary$sd_score
score_summary$upper_score <- score_summary$mean_score + score_summary$sd_score
smooth_score_series <- function(y_values) {
  smooth_days <- base::seq(base::min(days), base::max(days), length.out = 220)
  stats::spline(
    x = score_summary$day_numeric,
    y = y_values,
    xout = smooth_days,
    method = 'natural'
  )$y
}
score_curve <- base::data.frame(
  day_numeric = base::seq(base::min(days), base::max(days), length.out = 220),
  mean_score = smooth_score_series(score_summary$mean_score),
  lower_score = smooth_score_series(score_summary$lower_score),
  upper_score = smooth_score_series(score_summary$upper_score)
)
score_curve_lower <- base::pmin(score_curve$lower_score, score_curve$upper_score)
score_curve_upper <- base::pmax(score_curve$lower_score, score_curve$upper_score)
score_curve$lower_score <- score_curve_lower
score_curve$upper_score <- score_curve_upper

p_score_by_day <- ggplot2::ggplot(pca_df, ggplot2::aes(day_numeric, differentiation_score)) +
  ggplot2::geom_hline(
    yintercept = base::c(0, 1), linetype = "dashed",
    linewidth = gs(0.24), color = "black"
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
    linewidth = gs(0.52),
    color = 'black'
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = day_numeric),
    size = gs(0.95),
    alpha = 1
  ) +
  ggplot2::scale_x_continuous(breaks = days) +
  day_color_scale(name = 'Day') +
  ggplot2::labs(x = "Differentiation day", y = "Differentiation timing score") +
  theme_pub() +
  ggplot2::theme(
    legend.position = base::c(0.36, 0.78),
    legend.justification = base::c(0.5, 0.5),
    legend.direction = 'horizontal',
    legend.background = ggplot2::element_rect(fill = scales::alpha('white', 0.82), color = NA),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    legend.key.width = grid::unit(0.30, 'in'),
    legend.title = ggplot2::element_text(
      family = figure_family,
      face = 'bold',
      size = fs(6),
      hjust = 0.5
    )
  )

# 8.0 load cached leave-one-line-out validation -----------------

if (!file.exists(loo_validation_cache_path)) {
  stop(
    'Required leave-one-line-out validation cache is missing. Run ',
    'scripts/02_run_leave_one_line_out_validation.R before rendering.',
    call. = FALSE
  )
}

validation_source_paths <- c(
  'scripts/02_run_leave_one_line_out_validation.R',
  'functions/load_GSE122380_data.R',
  'functions/select_temporal_genes.R',
  'functions/score_differentiation_timing.R',
  'config/analysis.yml'
)
expected_validation_cache_key <- list(
  input_md5 = GSE122380_data$input_md5,
  implementation_md5 = unname(tools::md5sum(validation_source_paths)),
  package_versions = c(
    DESeq2 = as.character(utils::packageVersion('DESeq2')),
    edgeR = as.character(utils::packageVersion('edgeR'))
  ),
  expression_cpm_cutoff = expression_cpm_cutoff,
  lrt_padj_cutoff = lrt_padj_cutoff,
  vst_dynamic_range_cutoff = vst_dynamic_range_cutoff,
  n_pcs = n_timing_pcs,
  heldout_lines = sort(unique(as.character(meta$cell_line)))
)

message('CODEX_STEP loo: loading cached leave-one-line-out validation')
loo_cache <- readRDS(loo_validation_cache_path)
if (!isTRUE(identical(loo_cache$cache_key, expected_validation_cache_key))) {
  stop(
    'The leave-one-line-out validation cache is stale for the current inputs, ',
    'code, parameters, or package versions. Rerun ',
    'scripts/02_run_leave_one_line_out_validation.R.',
    call. = FALSE
  )
}

  loo_scores <- loo_cache$scores
  loo_polyline_all <- loo_scores[
    loo_scores$gene_set == 'All temporal',
    ,
    drop = FALSE
  ]
  if ('method' %in% base::names(loo_polyline_all)) {
    loo_polyline_all <- loo_polyline_all[
      loo_polyline_all$method == 'Polyline',
      ,
      drop = FALSE
    ]
  }
  loo_polyline_all$heldout_line <- base::factor(
    loo_polyline_all$heldout_line,
    levels = base::sort(base::unique(loo_polyline_all$heldout_line))
  )
  loo_display_excluded_lines <- '19190'
  loo_best_all <- loo_polyline_all[
    !loo_polyline_all$heldout_line %in% loo_display_excluded_lines,
    ,
    drop = FALSE
  ]
  loo_best_all$heldout_line <- base::factor(
    loo_best_all$heldout_line,
    levels = base::sort(base::unique(loo_best_all$heldout_line))
  )
  loo_line_levels <- base::levels(loo_best_all$heldout_line)
  loo_line_ncol <- 4L
  loo_line_nrow <- base::ceiling(base::length(loo_line_levels) / loo_line_ncol)
  loo_line_plots <- base::lapply(base::seq_along(loo_line_levels), function(i) {
    line_id <- loo_line_levels[[i]]
    line_df <- loo_best_all[loo_best_all$heldout_line == line_id, , drop = FALSE]
    row_index <- base::ceiling(i / loo_line_ncol)
    col_index <- ((i - 1L) %% loo_line_ncol) + 1L
    show_y_title <- col_index == 1L
    show_x_title <- row_index == loo_line_nrow
    line_r <- stats::cor(line_df$actual_day, line_df$predicted_day, use = 'complete.obs')
    line_r_label <- base::sprintf(
      "r = %.2f\nr\u00b2 = %.2f",
      line_r,
      line_r^2
    )

    ggplot2::ggplot(line_df, ggplot2::aes(actual_day, predicted_day)) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        color = 'black',
        linetype = 'dashed',
        linewidth = gs(0.24)
      ) +
      ggplot2::geom_point(
        color = 'black',
        size = gs(0.55),
        alpha = 1
      ) +
      ggplot2::geom_smooth(
        method = 'loess',
        formula = y ~ x,
        se = FALSE,
        span = 0.78,
        color = '#1260A4',
        linewidth = gs(0.52)
      ) +
      ggplot2::annotate(
        'text',
        x = base::min(days),
        y = base::max(days),
        label = line_r_label,
        hjust = 0,
        vjust = 1,
        family = figure_family,
        size = gfs(5.1),
        lineheight = 0.86
      ) +
      ggplot2::scale_x_continuous(breaks = days, limits = base::range(days)) +
      ggplot2::scale_y_continuous(breaks = days, limits = base::range(days)) +
      ggplot2::labs(
        title = base::paste('Line', line_id),
        x = if (show_x_title) 'Actual differentiation day' else NULL,
        y = if (show_y_title) 'Predicted differentiation day' else NULL
      ) +
      theme_pub() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = 'plain',
          hjust = 0,
          size = fs(7),
          family = figure_family,
          color = 'black'
        ),
        axis.text.x = ggplot2::element_text(size = fs(6), color = 'black', family = figure_family),
        axis.text.y = ggplot2::element_text(size = fs(6), color = 'black', family = figure_family),
        axis.title.x = ggplot2::element_text(size = fs(7), color = 'black', family = figure_family),
        axis.title.y = ggplot2::element_text(size = fs(7), color = 'black', family = figure_family),
        plot.margin = ggplot2::margin(5, 6, 6, 5)
      )
  })
  p_loo_line_predictions <- patchwork::wrap_plots(loo_line_plots, ncol = loo_line_ncol)

  loo_r <- stats::cor(loo_polyline_all$actual_day, loo_polyline_all$predicted_day, use = 'complete.obs')
  loo_residual <- loo_polyline_all$predicted_day - loo_polyline_all$actual_day
  loo_abs_error <- base::abs(loo_residual)
  loo_accuracy_summary <- base::data.frame(
    metric = base::c(
      'Held-out cell lines',
      'Held-out samples',
      'Correlation between actual and predicted day',
      'Squared Pearson correlation',
      'Mean absolute error',
      'Median absolute error',
      'Predictions within 1 day',
      'Predictions within 2 days'
    ),
    value = base::c(
      base::format(base::length(base::unique(loo_polyline_all$heldout_line)), big.mark = ','),
      base::format(base::nrow(loo_polyline_all), big.mark = ','),
      base::sprintf('%.3f', loo_r),
      base::sprintf('%.3f', loo_r^2),
      base::sprintf('%.2f days', base::mean(loo_abs_error, na.rm = TRUE)),
      base::sprintf('%.2f days', stats::median(loo_abs_error, na.rm = TRUE)),
      base::sprintf('%.1f%%', base::mean(loo_abs_error <= 1, na.rm = TRUE) * 100),
      base::sprintf('%.1f%%', base::mean(loo_abs_error <= 2, na.rm = TRUE) * 100)
    ),
    stringsAsFactors = FALSE
  )

  loo_lm_label <- format_lm_equation_label(
    x = loo_polyline_all$actual_day,
    y = loo_polyline_all$predicted_day
  )

  p_loo_predicted_vs_actual <- ggplot2::ggplot(
    loo_polyline_all,
    ggplot2::aes(actual_day, predicted_day)
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = 'black',
      linetype = 'dashed',
      linewidth = gs(0.24)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = heldout_line),
      size = gs(0.95),
      alpha = 1
    ) +
    ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = gs(0.45)
    ) +
    ggplot2::annotate(
      'text',
      x = base::min(loo_polyline_all$actual_day),
      y = base::max(loo_polyline_all$predicted_day),
      label = loo_lm_label,
      hjust = 0,
      vjust = 1,
      family = figure_family,
      size = gfs(5.5),
      lineheight = 0.95
    ) +
    ggplot2::scale_color_viridis_d(name = 'Held-out line', option = 'D') +
    ggplot2::scale_x_continuous(breaks = days) +
    ggplot2::scale_y_continuous(breaks = days) +
    ggplot2::labs(x = 'Actual differentiation day', y = 'Predicted differentiation day') +
    theme_pub() +
    ggplot2::theme(legend.position = 'none')

  p_loo_summary <- p_loo_predicted_vs_actual

  loo_timepoint_error <- loo_polyline_all
  loo_timepoint_error$absolute_error <- base::abs(loo_timepoint_error$residual)
  loo_timepoint_error$actual_day_factor <- base::factor(
    loo_timepoint_error$actual_day,
    levels = days
  )
  loo_timepoint_label <- stats::aggregate(
    absolute_error ~ actual_day_factor,
    data = loo_timepoint_error,
    FUN = function(error) {
      base::min(base::mean(error), stats::median(error))
    }
  )
  loo_timepoint_y <- stats::aggregate(
    absolute_error ~ actual_day_factor,
    data = loo_timepoint_error,
    FUN = function(error) {
      grDevices::boxplot.stats(error)$stats[5]
    }
  )
  names(loo_timepoint_label)[2] <- 'error_label'
  names(loo_timepoint_y)[2] <- 'label_y'
  loo_timepoint_label <- base::merge(
    loo_timepoint_label,
    loo_timepoint_y,
    by = 'actual_day_factor',
    sort = FALSE
  )
  loo_timepoint_label$label <- base::sprintf('%.1f', loo_timepoint_label$error_label)

  p_loo_timepoint_accuracy <- ggplot2::ggplot(
    loo_timepoint_error,
    ggplot2::aes(actual_day_factor, absolute_error)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = 'black',
      linewidth = gs(0.22),
      linetype = 'dashed'
    ) +
    ggplot2::geom_boxplot(
      width = 0.38,
      outlier.shape = NA,
      linewidth = gs(0.30),
      alpha = 0.70,
      staplewidth = 0.3,
      fill = '#BDBDBD',
      color = 'black'
    ) +
    ggplot2::geom_text(
      data = loo_timepoint_label,
      ggplot2::aes(actual_day_factor, label_y, label = label),
      inherit.aes = FALSE,
      vjust = -0.45,
      family = figure_family,
      size = gfs(4.7)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = base::c(0, 0.16))) +
    ggplot2::labs(
      x = 'Actual differentiation day',
      y = 'Absolute prediction error (days)'
    ) +
    ggplot2::theme_classic(base_size = fs(7), base_family = figure_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figure_family, color = 'black'),
      axis.title = ggplot2::element_text(size = fs(7)),
      axis.text.x = ggplot2::element_text(size = fs(6)),
      axis.text.y = ggplot2::element_text(size = fs(6)),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
message('CODEX_STEP figures: saving tutorial figures')

save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_reference_pca_and_day_correlation'),
  plot = p_reference_overview,
  width = 7.2,
  height = 3.81
)
save_drawn_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_temporal_heatmap'),
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
  width = 5.45,
  height = 3.35,
  png_scale = 2
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_temporal_clusters'),
  plot = p_cluster_trajectories,
  width = 7.2,
  height = 8.52
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_pca_day'),
  plot = p_pca_day,
  width = 7.2,
  height = 8.52
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_pc1_validation'),
  plot = p_pc1_validation,
  width = 7.2,
  height = 4.95
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_timing_polyline'),
  plot = p_timing_polyline,
  width = 7.2,
  height = 3.81
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_score_by_day'),
  plot = p_score_by_day,
  width = 7.2,
  height = 3.65
)

save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_loo_line_predictions'),
  plot = p_loo_line_predictions,
  width = 7.2,
  height = 5.75
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_loo_summary'),
  plot = p_loo_summary,
  width = 7.2,
  height = 3.15
)
save_figure_pair(
  path_stub = base::file.path(docs_figure_dir, 'GSE122380_loo_timepoint_accuracy'),
  plot = p_loo_timepoint_accuracy,
  width = 5.4,
  height = 2.75
)

base::message('CODEX_DONE built tutorial figures and summary tables')
