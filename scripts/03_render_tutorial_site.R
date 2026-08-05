# ----
# author:
# - Zoheb Khan
#
# script path:
# - scripts/03_render_tutorial_site.R
#
# input data:
# - tutorial/tutorial.Rmd
# - tutorial/_bookdown.yml
# - tutorial/_output.yml
# - tutorial/style.css
# - tutorial/assets/fonts/*.{otf,ttf}
# - scripts/01_build_tutorial_objects.R
# - tmp/GSE122380_leave_one_cell_line_out_validation.rds
#
# outputs:
# - docs/index.html
# - docs/style.css
# - docs/search_index.json
# - docs/reference-keys.txt
# - docs/libs/
# - docs/assets/figures/
# - docs/assets/fonts/
# - docs/.nojekyll
# ----

# 0.0 validate render inputs and dependencies -----------------

if (!requireNamespace('bookdown', quietly = TRUE)) {
  stop('The bookdown package is required to render the tutorial.', call. = FALSE)
}

source_font_dir = 'tutorial/assets/fonts'
site_fonts <- file.path(
  source_font_dir,
  c(
    'LatinModernSans-Bold.otf',
    'LatinModernSans-BoldOblique.otf',
    'LatinModernSans-Oblique.otf',
    'LatinModernSans-Regular.otf',
    'FiraCode-Retina.ttf'
  )
)
figure_fonts <- file.path(
  source_font_dir,
  c(
    'NimbusSans-Bold.otf',
    'NimbusSans-BoldItalic.otf',
    'NimbusSans-Italic.otf',
    'NimbusSans-Regular.otf'
  )
)

required_inputs <- c(
  'tutorial/tutorial.Rmd',
  'tutorial/_bookdown.yml',
  'tutorial/_output.yml',
  'tutorial/style.css',
  'scripts/01_build_tutorial_objects.R',
  'tmp/GSE122380_leave_one_cell_line_out_validation.rds',
  site_fonts,
  figure_fonts
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    'Missing required render inputs: ',
    paste(missing_inputs, collapse = ', '),
    '. Run scripts/02_run_leave_one_cell_line_out_validation.R before rendering.',
    call. = FALSE
  )
}

# 1.0 build tutorial objects and render the single source -----------------

tutorial_environment <- new.env(parent = baseenv())
source(
  'scripts/01_build_tutorial_objects.R',
  local = tutorial_environment
)

bookdown::render_book(
  input = 'tutorial',
  output_format = 'bookdown::gitbook',
  clean = TRUE,
  envir = tutorial_environment
)

required_site_outputs <- c(
  'docs/index.html',
  'docs/style.css',
  'docs/search_index.json',
  'docs/reference-keys.txt'
)
missing_site_outputs <- required_site_outputs[!file.exists(required_site_outputs)]
if (length(missing_site_outputs) > 0L) {
  stop(
    'Rendering did not create required site outputs: ',
    paste(missing_site_outputs, collapse = ', '),
    '.',
    call. = FALSE
  )
}

# 2.0 install static font assets and GitHub Pages marker -----------------

docs_font_dir <- 'docs/assets/fonts'
dir.create(docs_font_dir, recursive = TRUE, showWarnings = FALSE)
font_copy_ok <- file.copy(
  from = site_fonts,
  to = docs_font_dir,
  overwrite = TRUE,
  copy.date = TRUE
)
if (!all(font_copy_ok)) {
  stop('One or more tutorial font files could not be copied to docs.', call. = FALSE)
}

if (!file.create('docs/.nojekyll')) {
  stop('Could not create docs/.nojekyll.', call. = FALSE)
}

message('Rendered tutorial site: docs/index.html')
