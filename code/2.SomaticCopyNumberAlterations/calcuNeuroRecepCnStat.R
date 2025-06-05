
library('dplyr')


load(file = '/data/getPanCanNrCopyNumData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


nRCopyNumData <- rownames_to_column(nRCopyNumData, var = 'Approved.symbol') %>% 
  reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')

nRCopyNumData <- nRCopyNumData %>% left_join(tcgaPanCanSamples, by = join_by(SAMPLE_BARCODE)) %>% 
  left_join(neurotransmitterReceptors, by = join_by(Approved.symbol))


samsStat <- nRCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, cnStaus) %>% summarise(n = n()) %>% 
  group_by(DISEASE)  %>% summarize(numberOfPatient = n_distinct(SAMPLE_BARCODE))

##########
nRCnStatByDis <- nRCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, cnStaus) %>% count() %>% 
  group_by(DISEASE, cnStaus) %>% summarize(numberOfPatientAlt = n())
  
nRCnStatByDis <- nRCnStatByDis %>% subset(cnStaus != 0) %>% left_join(samsStat, by = join_by(DISEASE)) %>% 
  mutate(altFre = numberOfPatientAlt/numberOfPatient)

##########
nRCnStatByGroupDis <- nRCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, classesOfNeuroreceptors, cnStaus) %>% count() %>% 
  group_by(DISEASE, classesOfNeuroreceptors, cnStaus) %>% summarize(numberOfPatientAlt = n())

nRCnStatByGroupDis <- nRCnStatByGroupDis %>% subset(cnStaus != 0) %>% left_join(samsStat, by = join_by(DISEASE)) %>% 
  mutate(altFre = numberOfPatientAlt/numberOfPatient)

##########
nRCnStatByGeneDis <- nRCopyNumData %>% group_by(DISEASE, Approved.symbol, cnStaus) %>% summarize(numberOfPatientAlt = n())

nRCnStatByGeneDis <- nRCnStatByGeneDis %>% subset(cnStaus != 0) %>% left_join(samsStat, by = join_by(DISEASE)) %>% 
  mutate(altFre = numberOfPatientAlt/numberOfPatient)

##########

nRCnStatByGene <- nRCopyNumData %>% group_by(Approved.symbol, cnStaus) %>% summarize(numberOfPatientAlt = n()) %>% 
  subset(cnStaus != 0)

nRCnStatByGene$altFre <- nRCnStatByGene$numberOfPatientAlt/sum(samsStat$numberOfPatient)


save(nRCnStatByDis, nRCnStatByGroupDis, nRCnStatByGeneDis, nRCnStatByGene, file = '/data/neuroRecepCnStat.RData')


