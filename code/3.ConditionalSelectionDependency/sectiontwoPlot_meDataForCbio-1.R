library('dplyr')
library('tibble')

#######################
load(file = '/data/panCanNrCopyNumData.RData')

tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
tcgaPanCanCliData <- readRDS(file = '/data/tcgaPanCanCliData.rds')
tcgaPanCanCliData <- tcgaPanCanCliData %>% rownames_to_column(var = 'PATIENT_BARCODE')
tcgaPanCanSamples <- tcgaPanCanSamples %>% left_join(tcgaPanCanCliData, by ='PATIENT_BARCODE')
tcgaPanCanSamples <- tcgaPanCanSamples %>% dplyr::select(SAMPLE_BARCODE, DISEASE, age, gender, race, ajcc_stage, clinical_stage)

##########
nRCopyNumData <- rownames_to_column(nRCopyNumData, var = 'Approved.symbol') %>% 
  reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')

# nRCopyNumData <- rownames_to_column(proteinGeneCopyNumData[c('EGFR', 'TRIO'), ], var = 'Approved.symbol') %>%
#   reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')

# nRCopyNumData <- rownames_to_column(proteinGeneCopyNumData[c('EGFR', 'ARHGEF1'), ], var = 'Approved.symbol') %>%
#   reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')


nRCopyNumData <- nRCopyNumData %>% mutate(Type = 'CNA', cnStaus = ifelse(cnStaus == 2, 'AMP', ifelse(cnStaus == -2, 'HOMDEL', 'NEUTRAL'))) %>% 
  dplyr::select(SAMPLE_BARCODE, Approved.symbol, cnStaus, Type) %>% arrange(cnStaus)

##########
rhoCnStatByGene <- nRCopyNumData %>% subset(cnStaus != 'NEUTRAL') %>% group_by(Approved.symbol) %>% count() %>% arrange(desc(n))


showGenes <- rhoCnStatByGene$Approved.symbol[1:30]

# HTR3C,HTR3D,HTR3E,ADRB3,CHRNB3,CHRNA6,CHRNB2,ADRA1A,CHRNA2,CHRM3,OPRK1,OPRL1,HRH3,,CHRNA4,
# GABRA3,GABRQ,HTR2A,GABRE,GRIK2,GRIN2C,GLRA3,HRH4,GRM5,GRM7,GABRD,GABRR3,GRIN2B,GRIN3B,
# GRM6,DRD3


showGenesCopyNumData <- nRCopyNumData %>% mutate(cnStaus = ifelse(Approved.symbol %in% showGenes, cnStaus, 'NEUTRAL')) %>% arrange(cnStaus)

##########
setwd('/result/Section2/')

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

