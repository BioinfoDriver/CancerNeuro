
library('dplyr')
library('tibble')

##############
setwd('/data')
panCanMethyData <- data.table::fread(input = 'jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv')
panCanMethyData <- panCanMethyData %>% remove_rownames() %>% column_to_rownames(var = 'Composite Element REF')

# duplicated
panCanMethyData <- panCanMethyData[, !duplicated(substr(colnames(panCanMethyData), 1, 15))]
colnames(panCanMethyData) <- substr(colnames(panCanMethyData), 1, 15)


##############
tcgaPanCanSamples <- readRDS(file = 'tcgaPanCanSamples.rds')

turSmas <- intersect(colnames(panCanMethyData), tcgaPanCanSamples$SAMPLE_BARCODE)
normSams <- colnames(panCanMethyData)[substr(colnames(panCanMethyData), 14, 15) == 11]

pairdSams <- intersect(substr(turSmas, 1, 12), substr(normSams, 1, 12))
pairdTurSams <- turSmas[substr(turSmas, 1, 12) %in% pairdSams]
pairdNorSams <- normSams[substr(normSams, 1, 12) %in% pairdSams]


panCanTurMethy <- panCanMethyData %>% select(all_of(turSmas))
panCanPairdTurMethy <- panCanMethyData %>% select(all_of(pairdTurSams))
panCanPairdNormMethy <- panCanMethyData %>% select(all_of(pairdNorSams))


save(panCanTurMethy, panCanPairdTurMethy, panCanPairdNormMethy, file = 'panCanMethyData.RData')

