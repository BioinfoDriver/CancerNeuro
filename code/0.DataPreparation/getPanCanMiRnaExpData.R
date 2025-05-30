
library('dplyr')
library('tibble')
########################
setwd('/data')

panMiRnaExp <- read.csv(file = 'pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv', header = TRUE, sep = ',')
panMiRnaExp <- panMiRnaExp %>% subset(Correction == 'Corrected') %>% column_to_rownames(var = 'Genes') %>% mutate(Correction = NULL)
colnames(panMiRnaExp) <- substr(gsub(pattern = '.', replacement = '-', colnames(panMiRnaExp), fixed = TRUE), 1, 15)


# preprocess
panMiRnaExp[panMiRnaExp < 0] <- NA
panMiRnaExp <- log2(panMiRnaExp + 1)


#####
tcgaPanCanSamples <- readRDS(file = 'tcgaPanCanSamples.rds')

turSmas <- intersect(colnames(panMiRnaExp), tcgaPanCanSamples$SAMPLE_BARCODE)
normSams <- colnames(panMiRnaExp)[substr(colnames(panMiRnaExp), 14, 15) == 11]

pairdSams <- intersect(substr(turSmas, 1, 12), substr(normSams, 1, 12))
pairdTurSams <- turSmas[substr(turSmas, 1, 12) %in% pairdSams]
pairdNorSams <- normSams[substr(normSams, 1, 12) %in% pairdSams]


panCanTurMiRnaExp <- panMiRnaExp %>% select(all_of(turSmas))
panCanPairdTurMiRnaExp <- panMiRnaExp %>% select(all_of(pairdTurSams))
panCanPairdNormMiRnaExp <- panMiRnaExp %>% select(all_of(pairdNorSams))


save(panCanTurMiRnaExp, panCanPairdTurMiRnaExp, panCanPairdNormMiRnaExp, file = 'panCanMiRnaExpData.RData')

