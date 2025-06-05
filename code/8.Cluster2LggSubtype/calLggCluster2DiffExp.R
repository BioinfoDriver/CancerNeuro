library('dplyr')
library('tibble')

####################
load(file = '/data/panCanGeneExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


tcgaExpData <- panCanTurGeneExp
tcgaExpData <- tcgaExpData[as.character(neurotransmitterReceptors$NCBI.Gene.ID), ]
rownames(tcgaExpData) <- neurotransmitterReceptors$Approved.symbol  

####################
setwd('/result/section6/consensusCluster_New/New3')

lggCluster <- read.csv(file='panCanClusterNew3.k=5.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
lggCluster <- lggCluster %>% dplyr::rename(SAMPLE_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('PATIENT_BARCODE', 'SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE') %>% 
  arrange(Clusters, DISEASE) %>% mutate(Clusters = paste0('Cluster', Clusters)) %>% subset(DISEASE == 'LGG')


####################
lggExpData <- tcgaExpData[, intersect(colnames(tcgaExpData), lggCluster$SAMPLE_BARCODE)]
# Filtering
indices <- apply(lggExpData, 1, function(x) {(sum(is.na(x)) + sum(x == 0, na.rm = T))/length(x) < 0.2})
lggExpData <- lggExpData[indices, ]



lggExpData <- as.data.frame(t(lggExpData))
lggExpData <- lggExpData %>% rownames_to_column(var = "SAMPLE_BARCODE")
lggExpData <- merge(lggExpData, lggCluster[, c('SAMPLE_BARCODE', 'Clusters')], by = 'SAMPLE_BARCODE')


lggExpData <- lggExpData %>% column_to_rownames(var = 'SAMPLE_BARCODE') %>% mutate(Clusters = ifelse(Clusters== 'Cluster2', 'Poor', 'Good'))


difExpPvalue <- sapply(colnames(lggExpData)[-ncol(lggExpData)], function(geneID){
  
  panExp <- lggExpData[, c(geneID, 'Clusters')]
  colnames(panExp) <- c('geneID', 'Clusters')
  print(table(is.na(panExp$geneID)))
  p.value <- summary(aov(geneID ~ Clusters, panExp))[[1]]$'Pr(>F)'[1]
  
  return(p.value)
})

difExpQvalue <- p.adjust(difExpPvalue, method = 'fdr')


saveRDS(difExpQvalue, file = '/data/lggCluster2DiffExpGene.rds')
