

mc3MutData <- readRDS(file = '/data/mc3MutData.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')

neuroRecepMutData <- merge(mc3MutData, neurotransmitterReceptors, by.x = 'Hugo_Symbol', by.y = 'Approved.symbol')

neuroRecepMutData <- subset(neuroRecepMutData, 
  Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',  'In_Frame_Ins', 
   'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))



# Stat by patient
library('dplyr')
neuroRecepMutStatByDis <- neuroRecepMutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% 
    count(name = 'numOfMut')  %>% group_by(DISEASE)  %>% count(name = 'numOfPatWithMut')

neuroRecepMutStatByDisGroup <- neuroRecepMutData %>% group_by(DISEASE, classesOfNeuroreceptors, SAMPLE_BARCODE) %>% 
  count(name = 'numOfMut')  %>% group_by(DISEASE, classesOfNeuroreceptors)  %>% count(name = 'numOfPatWithMut')


disCount <- unique(neuroRecepMutData[, c('DISEASE', 'numOfPat')])

neuroRecepMutStatByDis <- merge(neuroRecepMutStatByDis, disCount, by = 'DISEASE')
neuroRecepMutStatByDisGroup <- merge(neuroRecepMutStatByDisGroup, disCount, by = 'DISEASE')

neuroRecepMutStatByDis <- neuroRecepMutStatByDis %>% mutate(altFre = numOfPatWithMut/numOfPat)
neuroRecepMutStatByDisGroup <- neuroRecepMutStatByDisGroup %>% mutate(altFre = numOfPatWithMut/numOfPat)



# Stat by gene
neuroRecepMutStatByDisGene <- neuroRecepMutData %>% group_by(DISEASE, SAMPLE_BARCODE, Hugo_Symbol) %>% count(name = 'numOfMut') %>% 
  group_by(DISEASE, Hugo_Symbol) %>% count(name = 'numOfPatWithMut') 


neuroRecepMutStatByDisGene <- merge(neuroRecepMutStatByDisGene, disCount, by = 'DISEASE')
neuroRecepMutStatByDisGene <- neuroRecepMutStatByDisGene %>% mutate(altFre = numOfPatWithMut/numOfPat)


neuroRecepMutStatByGene <- neuroRecepMutData %>% group_by(DISEASE, SAMPLE_BARCODE, Hugo_Symbol) %>% count(name = 'numOfMut') %>% 
  group_by(Hugo_Symbol) %>% count(name = 'numOfMut') 

neuroRecepMutStatByGene$altFre <- neuroRecepMutStatByGene$numOfMut/sum(disCount$numOfPat)


save(neuroRecepMutStatByDis, neuroRecepMutStatByDisGroup, neuroRecepMutStatByDisGene, neuroRecepMutStatByGene, file = '/data/neuroRecepMutStat.RData')



