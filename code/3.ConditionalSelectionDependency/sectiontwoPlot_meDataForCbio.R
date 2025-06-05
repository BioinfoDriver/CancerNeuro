library('dplyr')
library('tibble')

#######################
load(file = '/data/panCanNrCopyNumData.RData')
# copyNumData, nRCopyNumData


tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
tcgaPanCanCliData <- readRDS('/data/tcgaPanCanCliData.rds')
tcgaPanCanCliData <- tcgaPanCanCliData %>% rownames_to_column(var = 'PATIENT_BARCODE')
tcgaPanCanSamples <- tcgaPanCanSamples %>% left_join(tcgaPanCanCliData, by ='PATIENT_BARCODE')
tcgaPanCanSamples <- tcgaPanCanSamples %>% dplyr::select(SAMPLE_BARCODE, DISEASE, age, gender, race, ajcc_stage, clinical_stage)


##########
nRCopyNumData <- rownames_to_column(nRCopyNumData, var = 'Approved.symbol') %>% 
  reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')


# nRCopyNumData <- rownames_to_column(copyNumData[c('ERBB2', 'CHRM2'), ], var = 'Approved.symbol') %>%
#   reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')

# nRCopyNumData <- rownames_to_column(copyNumData[c('CCNE1', 'GRM5'), ], var = 'Approved.symbol') %>%
#   reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')



nRCopyNumData <- nRCopyNumData %>% mutate(Type = 'CNA', cnStaus = ifelse(cnStaus == 2, 'AMP', ifelse(cnStaus == -2, 'HOMDEL', 'NEUTRAL'))) %>% 
  dplyr::select(SAMPLE_BARCODE, Approved.symbol, cnStaus, Type) %>% arrange(cnStaus)

##########
rhoCnStatByGene <- nRCopyNumData %>% subset(cnStaus != 'NEUTRAL') %>% group_by(Approved.symbol) %>% count() %>% arrange(desc(n))

showGenes <- rhoCnStatByGene$Approved.symbol[1:24]
# showGenes <- rhoCnStatByGene %>% left_join(rhoFamilies, by = 'Approved.symbol') %>% group_by(Class) %>% slice_head(n = 6)
# showGenes <- showGenes$Approved.symbo

# ARHGAP39 ARHGAP30 DLC1 SRGAP2 ARHGAP33 ARHGAP4 ARHGDIA ARHGDIB ARHGDIG MCF2L2 ECT2 PLEKHG4B ARHGEF10
# TRIO ARHGEF2 RHOBTB2 RHOU RAC3 RHOD RAC1 RHOA


# showGenes <- c('ERBB2', 'CHRM2')
# showGenes <- c('CCNE1', 'GRM5')
showGenesCopyNumData <- nRCopyNumData %>% mutate(cnStaus = ifelse(Approved.symbol %in% showGenes, cnStaus, 'NEUTRAL')) %>% arrange(cnStaus)


##########
setwd('/result/section2/')

showGenesCopyNumData <- showGenesCopyNumData %>% mutate(Approved.symbol = ifelse(cnStaus == 'NEUTRAL', '', Approved.symbol), 
                                                        Type = ifelse(cnStaus == 'NEUTRAL', '', Type),
                                                        cnStaus = ifelse(cnStaus == 'NEUTRAL', '', cnStaus))

showGenesCopyNumData <- showGenesCopyNumData %>% subset(!(duplicated(SAMPLE_BARCODE) & Approved.symbol == ''))


tcgaPanCanSamples <- tcgaPanCanSamples %>% subset(SAMPLE_BARCODE %in% showGenesCopyNumData$SAMPLE_BARCODE) %>% 
  mutate(ajcc_stage = gsub(' ', '', ajcc_stage), clinical_stage = gsub(' ', '', clinical_stage))

tcgaPanCanSamples <- tcgaPanCanSamples %>% mutate_all(~ifelse(is.na(.), "N/A", .))


tcgaPanCanSamples <- tcgaPanCanSamples %>% rename(Cancer_Type = DISEASE, Age = age, Gender = gender, 
                                                  Race = race, Ajcc_Stage = ajcc_stage, Clinical_Stage = clinical_stage)

write.table(showGenesCopyNumData, file = 'cbioCnData.txt', sep = '\t', row.names = F, col.names = F, quote = F)
write.table(tcgaPanCanSamples, file = 'cbioSamAnno.txt', sep = '\t', row.names = F, col.names = T, quote = F)




#######################################ERBB2_CHRM2
showGenesCopyNumData <- showGenesCopyNumData %>% mutate(Approved.symbol = ifelse(cnStaus == 'HOMDEL', '', Approved.symbol)) %>% 
  mutate(Type = ifelse(cnStaus == 'HOMDEL', '', Type)) %>% mutate(cnStaus = ifelse(cnStaus == 'HOMDEL', '', cnStaus))


setwd('/result/section2/ERBB2_CHRM2')

write.table(showGenesCopyNumData, file = 'ERBB2_CHRM2_cbioCnData.txt', sep = '\t', row.names = F, col.names = F, quote = F)
write.table(tcgaPanCanSamples, file = 'ERBB2_CHRM2_cbioSamAnno.txt', sep = '\t', row.names = F, col.names = T, quote = F)


for(disease in unique(tcgaPanCanSamples$Cancer_Type)){
  
  disSams <- subset(tcgaPanCanSamples, Cancer_Type == disease)
  disCnData <- showGenesCopyNumData %>% mutate(SAMPLE_BARCODE = as.character(SAMPLE_BARCODE)) %>% 
    subset(SAMPLE_BARCODE %in% disSams$SAMPLE_BARCODE)
  
  table(disCnData$Approved.symbol, disCnData$cnStaus)
  write.table(disCnData, file = paste(disease, 'cbioCnData.txt', sep = '_'), sep = '\t', row.names = F, col.names = F, quote = F)
  write.table(disSams, file = paste(disease, 'cbioSamAnno.txt', sep = '_'), sep = '\t', row.names = F, col.names = T, quote = F)
  
}
table(merge(tcgaPanCanSamples, showGenesCopyNumData, by = 'SAMPLE_BARCODE')[, c('Approved.symbol', 'cnStaus', 'Cancer_Type')])


#######################################CCNE1_GRM5
showGenesCopyNumData <- showGenesCopyNumData %>% mutate(Approved.symbol = ifelse(cnStaus == 'HOMDEL', '', Approved.symbol)) %>% 
  mutate(Type = ifelse(cnStaus == 'HOMDEL', '', Type)) %>% mutate(cnStaus = ifelse(cnStaus == 'HOMDEL', '', cnStaus))


setwd('/result/section2/CCNE1_GRM5')

write.table(showGenesCopyNumData, file = 'CCNE1_GRM5_cbioCnData.txt', sep = '\t', row.names = F, col.names = F, quote = F)
write.table(tcgaPanCanSamples, file = 'CCNE1_GRM5_cbioSamAnno.txt', sep = '\t', row.names = F, col.names = T, quote = F)


for(disease in unique(tcgaPanCanSamples$Cancer_Type)){
  
  disSams <- subset(tcgaPanCanSamples, Cancer_Type == disease)
  disCnData <- showGenesCopyNumData %>% mutate(SAMPLE_BARCODE = as.character(SAMPLE_BARCODE)) %>% 
    subset(SAMPLE_BARCODE %in% disSams$SAMPLE_BARCODE)
  table(disCnData$Approved.symbol, disCnData$cnStaus)
  
  write.table(disCnData, file = paste(disease, 'cbioCnData.txt', sep = '_'), sep = '\t', row.names = F, col.names = F, quote = F)
  write.table(disSams, file = paste(disease, 'cbioSamAnno.txt', sep = '_'), sep = '\t', row.names = F, col.names = T, quote = F)
  
}
table(merge(tcgaPanCanSamples, showGenesCopyNumData, by = 'SAMPLE_BARCODE')[, c('Approved.symbol', 'cnStaus', 'Cancer_Type')])
