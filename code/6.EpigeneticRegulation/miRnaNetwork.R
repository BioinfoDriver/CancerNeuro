
library(dplyr)

miRnaTarByDisease <- readRDS(file = 'D:/CancerNeuroscience/Github/data/nROfmiRNARegulatorFLM.rds')

miRnaTarByDisease <- subset(miRnaTarByDisease, pValue < 0.05 & Estimate < -0.15)

network <- miRnaTarByDisease %>% group_by(miRNA, nRgene) %>% summarise(width = n()) %>% ungroup() %>% arrange(desc(width))

miRnaDegree <- miRnaTarByDisease %>% select(miRNA, nRgene) %>% distinct(miRNA, nRgene) %>% 
  group_by(miRNA) %>% count(name = 'degree') %>% mutate(nodeType = 1) %>% dplyr::rename(node = miRNA) %>% arrange(desc(degree))

targetDegree <- miRnaTarByDisease %>% select(miRNA, nRgene) %>% distinct(miRNA, nRgene) %>% 
  group_by(nRgene) %>% count(name = 'degree') %>% mutate(nodeType = 2) %>% dplyr::rename(node = nRgene)

nodeAttributes <- rbind.data.frame(miRnaDegree, targetDegree)

networkFilter <- subset(network, width > 2)


setwd('D:/CancerNeuroscience/Github/result/section4')
write.table(network, file = 'miRnaTargerNetworkTargetScan.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(nodeAttributes, file = 'networkNodeAttributesTargetScan.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(networkFilter, file = 'miRnaTargerNetworkTargetScanFilter.txt', sep = '\t', col.names = T, row.names = F, quote = F)

