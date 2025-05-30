
geneInfo <- readRDS(file = '/data/geneInfo.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


neurotransmitterReceptors <- merge(neurotransmitterReceptors, 
                                   select(geneInfo, gene_name, exon_breaks), by.x = 'Approved.symbol', by.y = 'gene_name')

neurotransmitterReceptors$exon_breaks <- as.character(neurotransmitterReceptors$exon_breaks)
geneInfo$exon_breaks <- as.character(geneInfo$exon_breaks)

geneInfo <- subset(geneInfo, exon_breaks %in% neurotransmitterReceptors$exon_breaks)



# random gene set (all)
geneGroup <- table(neurotransmitterReceptors$exon_breaks)

set.seed(1314)
randomGenesList <- lapply(1:1000, function(i){
  
  randomGenes <- lapply(seq(length(geneGroup)), function(j){
    
    gGroup <- geneGroup[j]
    gInfo <- subset(geneInfo, exon_breaks %in% names(gGroup))
    gInfo <- gInfo[sample(nrow(gInfo), size = gGroup), c('gene_name', 'exon_breaks')]
    
    return(gInfo)
  })
  
  randomGenes <- do.call(rbind, randomGenes)
  
  return(randomGenes)
})



# random gene set (by group)
geneGroup <- table(neurotransmitterReceptors$classesOfNeuroreceptors, neurotransmitterReceptors$exon_breaks)

set.seed(1314)
randomGenesListByGroup <- lapply(seq(nrow(geneGroup)), function(i){
  
  gGroup <- geneGroup[i, ]
  gGroup <- gGroup[gGroup > 0]
  
  randomGenesByGroup <- lapply(1:1000, function(j){
    
    randomGenes <- lapply(seq(length(gGroup)), function(k){
      
      ggGroup <- gGroup[k]
      
      gInfo <- subset(geneInfo, exon_breaks %in% names(ggGroup))
      gInfo <- gInfo[sample(nrow(gInfo), size = ggGroup), c('gene_name', 'exon_breaks')]
      
      gInfo$classesOfNeuroreceptors <- rownames(geneGroup)[i]
      return(gInfo)
    })
    
    randomGenes <- do.call(rbind, randomGenes)
    randomGenes <- as.data.frame(randomGenes)
  })
  
  return(randomGenesByGroup)
})


save(randomGenesList, randomGenesListByGroup, file = '/data/randomGenesListBasedOnExonLength.RData')

