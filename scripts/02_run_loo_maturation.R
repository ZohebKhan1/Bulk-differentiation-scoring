library(DESeq2)
library(edgeR)

source('data/load_GSE122380_data.R')
source('functions/run_loo_maturation.R')
source('functions/get_temporal_genes.R')
source('functions/score_differentiation_timing.R')

validation_cache_path <- 'cache/GSE122380_loo_maturation.rds'

GSE122380_data <- load_GSE122380_data()
validation_output <- run_loo_maturation(
  counts = GSE122380_data$counts,
  vst = GSE122380_data$vst,
  metadata = GSE122380_data$metadata
)

dir.create(dirname(validation_cache_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(validation_output, validation_cache_path)

validation_output$summary
message('Saved optional validation cache: ', validation_cache_path)
