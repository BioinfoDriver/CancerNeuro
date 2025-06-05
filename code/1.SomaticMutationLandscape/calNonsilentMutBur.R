library('dplyr')
###################
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

mc3MutData <- subset(mc3MutData, Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',  'In_Frame_Ins', 
                     'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))


meanOfNonsilentMutBurden <- mc3MutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% count(name = 'numOfNonsiletMut') %>% 
  group_by(DISEASE) %>% summarise(n = mean(numOfNonsiletMut))


saveRDS(meanOfNonsilentMutBurden, file = '/data/meanOfNonsilentMutBurden.rds')
