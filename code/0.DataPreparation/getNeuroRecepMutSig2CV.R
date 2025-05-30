

neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


setwd('/data/Bailey_Cell_2018_CancerDriverGenes/MutSig2CV')

list.files <- list.files()
list.files <- list.files[!grepl("^(README|UCEChyb)\\.txt$", list.files)]

mutSigGeneQvalue <- lapply(list.files, function(file.name){
  
  mutSigRes <- read.csv(file = file.name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  mutSigRes <- subset(mutSigRes, gene %in% neurotransmitterReceptors$Approved.symbol)
  
  mutSigRes <- mutSigRes[, c('gene', 'qvalue')]
  colnames(mutSigRes) <- c('gene', paste0(gsub('.txt', '', file.name), c('qvalue')))
  
  
  return(mutSigRes)
})


mutSigGeneQvalue <- Reduce(x = mutSigGeneQvalue, function(x, y) merge(x, y, by = 'gene'))


mutSigGenePvalue <- lapply(list.files, function(file.name){
  
  mutSigRes <- read.csv(file = file.name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  mutSigRes <- subset(mutSigRes, gene %in% neurotransmitterReceptors$Approved.symbol)
  
  mutSigRes <- mutSigRes[, c('gene', 'pvalue')]
  colnames(mutSigRes) <- c('gene', paste0(gsub('.txt', '', file.name), c('pvalue')))
  
  
  return(mutSigRes)
})


mutSigGenePvalue <- Reduce(x = mutSigGenePvalue, function(x, y) merge(x, y, by = 'gene'))


save(mutSigGenePvalue, mutSigGeneQvalue, file = '/data/neuroRecepMutSig2CV.RData')

