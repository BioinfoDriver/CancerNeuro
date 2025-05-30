
library('data.table')
library('dplyr')


setwd('/data')
tcgaPanCanSamples <- readRDS(file = 'tcgaPanCanSamples.rds')

panCanGeneExp <- fread(input = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv', sep = '\t')


geneInfo <- as.data.frame(do.call(rbind, strsplit(panCanGeneExp$gene_id, "\\|")))
colnames(geneInfo) <- c('geneName', 'geneID')

panCanGeneExp <- panCanGeneExp %>% as.data.frame()
rownames(panCanGeneExp) <- geneInfo$geneID


panCanGeneExp <- panCanGeneExp %>% mutate(gene_id  = NULL)
panCanGeneExp <- panCanGeneExp %>% select(-which(duplicated(substr(colnames(panCanGeneExp), 1, 15))))
colnames(panCanGeneExp) <- substr(colnames(panCanGeneExp), 1, 15)


turSmas <- intersect(colnames(panCanGeneExp), tcgaPanCanSamples$SAMPLE_BARCODE)
normSams <- colnames(panCanGeneExp)[substr(colnames(panCanGeneExp), 14, 15) == 11]

pairdSams <- intersect(substr(turSmas, 1, 12), substr(normSams, 1, 12))
pairdTurSams <- turSmas[substr(turSmas, 1, 12) %in% pairdSams]
pairdNorSams <- normSams[substr(normSams, 1, 12) %in% pairdSams]



# preprocess
panCanGeneExp[panCanGeneExp < 0] <- NA
panCanGeneExp <- log2(panCanGeneExp + 1)


panCanTurGeneExp <- panCanGeneExp %>% select(all_of(turSmas))
panCanPairdTurGeneExp <- panCanGeneExp %>% select(all_of(pairdTurSams))
panCanPairdNormGeneExp <- panCanGeneExp %>% select(all_of(pairdNorSams))


save(geneInfo, panCanTurGeneExp, panCanPairdTurGeneExp, panCanPairdNormGeneExp, file = 'panCanGeneExpData.RData')



