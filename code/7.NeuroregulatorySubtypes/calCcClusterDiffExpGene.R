library('dplyr')
library('tibble')

####################
load(file = '/data/panCanGeneExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


setwd('/result/result6/consensusCluster_New/New3')

panCanCluster <- read.csv(file='panCanClusterNew3.k=5.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% rename(SAMPLE_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE') %>% arrange(Clusters, DISEASE) %>% 
  mutate(Clusters = paste0('Cluster', Clusters))


panCanTurGeneExp <- as.data.frame(t(panCanTurGeneExp))
panCanTurGeneExp <- panCanTurGeneExp %>% rownames_to_column(var = "SAMPLE_BARCODE")
panCanTurGeneExp <- merge(panCanTurGeneExp, panCanCluster[, c('SAMPLE_BARCODE', 'Clusters')], by = 'SAMPLE_BARCODE')


panCanTurGeneExp <- panCanTurGeneExp[, intersect(colnames(panCanTurGeneExp), 
                                                 c('SAMPLE_BARCODE', 'Clusters', as.character(neurotransmitterReceptors$NCBI.Gene.ID)))] %>% 
  column_to_rownames(var = 'SAMPLE_BARCODE')



difExpPvalue <- sapply(colnames(panCanTurGeneExp)[-112], function(geneID){
  
  panExp <- panCanTurGeneExp[, c(geneID, 'Clusters')]
  colnames(panExp) <- c('geneID', 'Clusters')
  p.value <- summary(aov(geneID ~ Clusters, panExp))[[1]]$'Pr(>F)'[1]
  
  return(p.value)
})

difExpQvalue <- p.adjust(difExpPvalue, method = 'fdr')


names(difExpQvalue) <- neurotransmitterReceptors$Approved.symbol[match(names(difExpQvalue), neurotransmitterReceptors$NCBI.Gene.ID)]

saveRDS(difExpQvalue, file = '/data/ccClusterDiffExpGene_1.rds')

# 
# maxExpEachGroup <- panCanTurGeneExp %>% group_by(Clusters) %>% summarise(across(.fns = ~ mean(., na.rm = TRUE))) %>% as.data.frame()
# colnames(maxExpEachGroup)[-1] <- neurotransmitterReceptors$Approved.symbol[match(names(maxExpEachGroup)[-1], neurotransmitterReceptors$NCBI.Gene.ID)]
# 
# clusterMaxExpGene <- maxExpEachGroup %>% mutate(Clusters = NULL) %>% summarise(across(.fns = ~which.max(.))) %>% t() %>% as.data.frame %>%
#   merge.data.frame(neurotransmitterReceptors, by = "row.names") %>% arrange(V1) %>% select(Approved.symbol, V1, Chromosome)
# 
# saveRDS(clusterMaxExpGene, file = '/data/clusterMaxExpGene.rds')
# 

