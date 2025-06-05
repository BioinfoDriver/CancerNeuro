
library(dplyr)
library(reshape2)
#########################

geneInfo <- readRDS(file = '/data/geneInfo.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')

neurotransmitterReceptors <- merge(neurotransmitterReceptors, 
                                   select(geneInfo, gene_name, full_length, breaks), by.x = 'Approved.symbol', by.y = 'gene_name')


write.table(neurotransmitterReceptors, file = '/result/section1/nrReceptors.txt',
          sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)


##########
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

disCount <- mc3MutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% count(name = 'samMutNum') %>% 
  group_by(DISEASE) %>% count(name = 'numOfPat')


write.table(disCount, file = '/result/section1/numberOfPatWithMutData.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

##########
gainLoFMutStat <- readRDS(file = '/data/gainLoFMutStat.rds')

write.table(gainLoFMutStat, file = '/result/section1/gainLoFMutStat.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)


##########
load(file = '/data/neuroRecepMutStat.RData')
load(file = '/data/neuroRecepMutSig2CV.RData')


mutSigGenePvalue <- melt(mutSigGenePvalue, id.vars = 'gene', value.name = 'pValue', variable.name = 'cancerType')
mutSigGeneQvalue <- melt(mutSigGeneQvalue, id.vars = 'gene', value.name = 'qValue', variable.name = 'cancerType')

mutSigGenePvalue <- mutSigGenePvalue %>% mutate(cancerType = gsub('pvalue', '', cancerType))
mutSigGeneQvalue <- mutSigGeneQvalue %>% mutate(cancerType = gsub('qvalue', '', cancerType))


# mutSigGenePvalue <- subset(mutSigGenePvalue, pValue < 0.01)

sigMut <- merge(mutSigGenePvalue, mutSigGeneQvalue, by = c('gene', 'cancerType'))


neuroRecepMutStatByGene$DISEASE = 'PANCAN'
neuroRecepMutStatByGene$numOfPat = 8217
neuroRecepMutStatByGene <- rename(neuroRecepMutStatByGene, numOfPatWithMut = numOfMut)
neuroRecepMutStatByGene <- neuroRecepMutStatByGene[, colnames(neuroRecepMutStatByDisGene)]

nrMutSat <- rbind.data.frame(neuroRecepMutStatByGene, neuroRecepMutStatByDisGene)

sigMut <- rename(sigMut, Hugo_Symbol = gene, DISEASE = cancerType)
sigMut <- merge(sigMut, nrMutSat, by = c('Hugo_Symbol', 'DISEASE'))


write.table(sigMut, file = '/result/section1/nrReceptorsMutFreq.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)

##########
load(file = 'E:/CancerNeuroscience/data/neuroRecepMutComToRandom.RData')

comResByDisGroupPvalue <- comResByDisGroupPvalue %>% as.data.frame() %>% rownames_to_column(var = 'DISEASE')

comResByDisGroupPvalue <- melt(comResByDisGroupPvalue, id.vars = 'DISEASE', 
                               variable.name = 'classesOfNeuroreceptors', value.name = 'pValue')


neuroRecepMutStatByDisGroup <- merge(neuroRecepMutStatByDisGroup, comResByDisGroupPvalue, 
                                     by = c('DISEASE', 'classesOfNeuroreceptors'))


write.table(neuroRecepMutStatByDisGroup, file = '/result/section1/nrMutStatByDisGroup.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)



##########
load('/data/neuroRecepMutComToRandomExonLength_Threshold_LargeThan.RData')
load(file = '/data/neuroRecepMutStat.RData')

comResByDisGroupPvalue <- comResByDisGroupPvalue %>% as.data.frame() %>% rownames_to_column(var = 'DISEASE')

comResByDisGroupPvalue <- melt(comResByDisGroupPvalue, id.vars = 'DISEASE', 
                               variable.name = 'classesOfNeuroreceptors', value.name = 'pValue')


neuroRecepMutStatByDisGroup <- merge(neuroRecepMutStatByDisGroup, comResByDisGroupPvalue, 
                                     by = c('DISEASE', 'classesOfNeuroreceptors'))


write.table(neuroRecepMutStatByDisGroup, file = '/result/section1/nrMutStatByDisGroupExonLength_Threshold_LargeThan.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)



