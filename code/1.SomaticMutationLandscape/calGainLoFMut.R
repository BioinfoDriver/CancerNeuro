
library(dplyr)

#############
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


neuroRecepMutData <- merge(mc3MutData, neurotransmitterReceptors, by.x = 'Hugo_Symbol', by.y = 'Approved.symbol')

neuroRecepMutData <- subset(neuroRecepMutData, 
                            Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
                                                          'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))

neuroRecepMutDataByGene <- plyr::dlply(neuroRecepMutData, .variables = 'Hugo_Symbol')


# Across the pan-cancer cohort, hotspot mutations were defined as missense or in-frame mutations
# at the same protein amino acid in > 2 patient samples.
# The fraction of hotspot mutations per gene was calculated as the total number of hotspot mutations
# over the total number of non-silent mutations found in that gene.

# The fraction of LoF mutations (defined as Frame_Shift_Ins, Frame_Shift_Del, Nonsense_Mutation, Nonstop_Mutation,
# Splice_Site, and Tanslation_Start_Site) per gene was calculated as the total number of LoF mutations over the total number of
# non-silent mutations in that gene.
#
# Genes with > 30% hotspot mutations, < 20% LoF mutations, and ≥ 5 unique hotspot mutation
# positions were identified as enriched with hotspot mutations, while genes containing > 30% LoF mutations, < 30% hotspot
# mutations, and ≥ 10 LoF mutations were identified as enriched with LoF mutations.

fracOfLoF <- lapply(neuroRecepMutDataByGene, function(geneMutData){
  
  lofMut <- subset(geneMutData, Variant_Classification %in% c('Frame_Shift_Ins', 'Frame_Shift_Del', 'Nonsense_Mutation', 'Nonstop_Mutation',
                                                              'Splice_Site', 'Translation_Start_Site'))
  
  numOfLofMut <- nrow(lofMut)
  numOfNonsMut <- nrow(geneMutData)
  lofMf <- numOfLofMut/numOfNonsMut
  
  res <- c(numOfLofMut, lofMf, numOfNonsMut)
  
  names(res) <- c('numOfLofMut', 'lofMf', 'numOfNonsMut')
  return(res)
})

fracOfLoF <- do.call(rbind, fracOfLoF)


fracOfGainMut <- lapply(neuroRecepMutDataByGene, function(geneMutData){
  
  gainMut <- subset(geneMutData, Variant_Classification %in% c('In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation'))
  
  gainMut$mutPos <- gainMut$HGVSp
  gainMut$mutPos <- sapply(strsplit(gainMut$mutPos, split = '_'), function(x) x[1])
  gainMut$mutPos <- gsub("[A-Za-z]+$", '', gainMut$mutPos)
  
  
  mutPosi <- table(gainMut$mutPos)
  numOfgainMut <- sum(mutPosi[mutPosi >= 2])
  numOfHotMutPosi <- length(mutPosi[mutPosi >= 2])
  
  numOfNonsMut <- nrow(geneMutData)
  
  gainMf <- numOfgainMut/numOfNonsMut
  
  res <- c(numOfgainMut, numOfHotMutPosi, gainMf, numOfNonsMut)
  
  names(res) <- c('numOfgainMut', 'numOfHotMutPosi', 'gainMf', 'numOfNonsMut')
  return(res)
})

fracOfGainMut <- do.call(rbind, fracOfGainMut)


gainLoFMutStat <- cbind.data.frame(fracOfGainMut[, c('numOfgainMut', 'numOfHotMutPosi', 'gainMf')], fracOfLoF)
saveRDS(gainLoFMutStat, file = '/data/gainLoFMutStat.rds')

