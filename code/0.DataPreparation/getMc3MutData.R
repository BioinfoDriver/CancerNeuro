library('data.table')

setwd('/data/')
tcgaPanCanSamples <- readRDS(file = 'tcgaPanCanSamples.rds')

mc3MutData <- fread(file = 'mc3.v0.2.8.PUBLIC.maf', stringsAsFactors = FALSE)
mc3MutData$SAMPLE_BARCODE <- substr(mc3MutData$Tumor_Sample_Barcode, 1, 15)


mc3MutData <- merge(x = mc3MutData, y = tcgaPanCanSamples, by = 'SAMPLE_BARCODE')
mc3MutData <- subset(mc3MutData, FILTER == 'PASS' | (DISEASE %in% c('LAML', 'OV') & FILTER == 'wga'))


mutCount <- dplyr::count(mc3MutData, SAMPLE_BARCODE, name = 'samMutNum')
mutCount <- subset(mutCount, samMutNum <= 1000)

mc3MutData <- merge(x = mc3MutData, y = mutCount, by = 'SAMPLE_BARCODE')

mutCount <- dplyr::count(mc3MutData, DISEASE, SAMPLE_BARCODE, name = 'samMutNum')
disCount <- dplyr::count(mutCount, DISEASE, name = 'numOfPat')

mc3MutData <- merge(x = mc3MutData, y = disCount, by = 'DISEASE')

saveRDS(mc3MutData, file = 'mc3MutData.rds')
