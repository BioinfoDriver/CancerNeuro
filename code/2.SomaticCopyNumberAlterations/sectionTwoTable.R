
library('data.table')
library('dplyr')

########
load(file = '/data/neuroRecepCnStat.RData')
load(file = '/data/neuroRecepCnComToRandom.RData')


########
comResByDisGroupPvalue <- comResByDisGroupPvalue %>% melt() %>%  
  rename(DISEASE = Var1, classesOfNeuroreceptors = Var2, pValue = value) %>% mutate(DISEASE = as.character(DISEASE)) 

comResByDisGroupPvalue$cnStaus <- as.numeric(do.call(rbind, strsplit(comResByDisGroupPvalue$DISEASE, split = '_'))[, 2])
comResByDisGroupPvalue$DISEASE <- do.call(rbind, strsplit(comResByDisGroupPvalue$DISEASE, split = '_'))[, 1]

nRCnStatByGroupDis <- nRCnStatByGroupDis %>% right_join(comResByDisGroupPvalue, by = c('DISEASE', 'classesOfNeuroreceptors', 'cnStaus')) 
nRCnStatByGroupDis <- nRCnStatByGroupDis %>% mutate(cnStaus = ifelse(cnStaus == 2, 'Gain', 'Loss'))


write.table(nRCnStatByGroupDis, file = '/result/section2/nRCnStatByGroupDis.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)

########
comResByDisPvalue <- comResByDisPvalue %>% select(DISEASE, cnStaus, numberOfPatientAlt, numberOfPatient, altFre, pValue) %>% 
  mutate(cnStaus = ifelse(cnStaus == 2, 'Gain', 'Loss')) %>% arrange(DISEASE, cnStaus)

write.table(comResByDisPvalue, file = '/result/section2/nRCnStatByDis.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)

########

nRCnStatByGene <- nRCnStatByGene %>% mutate(cnStaus = ifelse(cnStaus == 2, 'Gain', 'Loss'))

write.table(nRCnStatByGene, file = '/result/section2/nRCnStatByGene.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)

########
nRCnStatByGeneDis <- nRCnStatByGeneDis %>% mutate(cnStaus = ifelse(cnStaus == 2, 'Gain', 'Loss'))

write.table(nRCnStatByGeneDis, file = '/result/section2/nRCnStatByGeneDis.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)


########
nrGisticPeakQvalue <- readRDS(file = '/data/nrGisticPeakQvalue.rds')

nrGisticPeakQvalue <- nrGisticPeakQvalue %>% as.data.frame() %>% 
  subset(!(Disease %in% c('COADREAD', 'GBMLGG', 'KIPAN', 'STES'))) %>% 
  select(Approved.symbol, Unique.Name:Wide.Peak.Limits, Direction, q.values, Disease, Direction)


write.table(nrGisticPeakQvalue, file = '/result/section2/nrGisticPeakQvalue.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)

########
enrichResideInPeak <- readRDS(file = '/data/enrichResideInPeak.rds')

write.table(enrichResideInPeak, file = '/result/section2/enrichResideInPeak.txt', sep = '\t', col.names = T, row.names = F, quote = F)
