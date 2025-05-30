
library('readxl')

setwd('/data/SanchezVega_Cell_2018_OncogenicSignalingPathways')

tcgaPanCanSamples <- read_xlsx(path = 'mmc1.xlsx', skip = 2, col_names = TRUE)
saveRDS(tcgaPanCanSamples, file = '/data/tcgaPanCanSamples.rds')