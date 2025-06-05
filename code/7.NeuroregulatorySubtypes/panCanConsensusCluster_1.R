
library('tibble')
library('ConsensusClusterPlus')
library('dplyr')


load(file = '/data/panCanGeneExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


# Filtering
panCanTurGeneExp <- panCanTurGeneExp[as.character(neurotransmitterReceptors$NCBI.Gene.ID), ]
rownames(panCanTurGeneExp) <- neurotransmitterReceptors$Approved.symbol

# anaGenes <- names(sort(apply(panCanTurGeneExp, 1, sd, na.rm = T), decreasing = T))[1:56]
# panCanTurGeneExp <- panCanTurGeneExp[anaGenes, ]


# Z-normalization
panCanTurGeneExp <- as.data.frame(t(panCanTurGeneExp))
panCanTurGeneExp <- panCanTurGeneExp %>% rownames_to_column(var = "SAMPLE_BARCODE")
panCanTurGeneExp <- merge(panCanTurGeneExp, tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE')


panCanTurGeneExp <- split.data.frame(panCanTurGeneExp, f = panCanTurGeneExp$DISEASE)

panCanTurGeneExp <- lapply(panCanTurGeneExp, function(geneExp){
  
  geneExp <- geneExp %>% mutate(DISEASE = NULL) %>% remove_rownames() %>% column_to_rownames(var = 'SAMPLE_BARCODE')
  
  geneExp <- t(scale(geneExp, center = TRUE, scale = TRUE))
  return(geneExp)
})

panCanTurGeneExp <- do.call(cbind, panCanTurGeneExp)



# Consensus Cluster
dt = as.dist(1-cor(panCanTurGeneExp, method="pearson", use = "pairwise.complete.obs"))

setwd('/results/section6/')

# ccRes <- ConsensusClusterPlus(
  # d=dt, maxK = 9, reps=100, pItem=0.8, pFeature=1, clusterAlg="hc", 
  # title="panCanClusterNew2",
  # innerLinkage="average", finalLinkage="ward.D2", distance="pearson", ml=NULL,
  # tmyPal=NULL, seed=1024, plot='png', writeTable=TRUE,
  # weightsItem=NULL, weightsFeature=NULL, verbose=TRUE, corUse="pairwise.complete.obs")

# saveRDS(ccRes, file = '/result/section6/panCanClusterNew2/ccRes.rds')



ccRes <- ConsensusClusterPlus(
  d=dt, maxK = 9, reps=100, pItem=0.8, pFeature=1, clusterAlg="hc", 
  title="panCanClusterNew3",
  innerLinkage="average", finalLinkage="average", distance="pearson", ml=NULL,
  tmyPal=NULL, seed=1024, plot='png', writeTable=TRUE,
  weightsItem=NULL, weightsFeature=NULL, verbose=TRUE, corUse="pairwise.complete.obs")

saveRDS(ccRes, file = '/result/section6/panCanClusterNew3/ccRes.rds')

