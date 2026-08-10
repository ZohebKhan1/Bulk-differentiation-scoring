library(bookdown)

source_font_dir = 'docs/tutorial/assets/fonts'
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
# 1.0 build tutorial objects and render the single source -----------------

source('scripts/01_build_tutorial_objects.R')

render_book(
  input = 'docs/tutorial',
  output_format = 'bookdown::gitbook',
  clean = TRUE,
  envir = globalenv()
)
unlink('docs/reference-keys.txt')

# 2.0 install static font assets and GitHub Pages marker -----------------

docs_font_dir <- 'docs/assets/fonts'
dir.create(docs_font_dir, recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(
  from = site_fonts,
  to = docs_font_dir,
  overwrite = TRUE,
  copy.date = TRUE
))
invisible(file.create('docs/.nojekyll'))

message('Rendered tutorial site: docs/index.html')
