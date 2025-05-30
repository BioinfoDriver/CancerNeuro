
setwd('/data/NeuroReceptors')

list.files <- list.files(pattern = 'Receptors.txt')

neurotransmitterReceptors <- sapply(list.files, function(file.name){
  
  neuroreceptors <- read.csv(file = file.name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  neuroreceptors$classesOfNeuroreceptors <- gsub(pattern = 's.txt', replacement = '', x = file.name)
  
  return(neuroreceptors)
}, simplify = FALSE)


neurotransmitterReceptors <- do.call(rbind, neurotransmitterReceptors)
rownames(neurotransmitterReceptors) <- neurotransmitterReceptors$Approved.symbol

neurotransmitterReceptors <- subset(neurotransmitterReceptors, Locus.type != 'pseudogene')

saveRDS(neurotransmitterReceptors, file = '/data/neurotransmitterReceptors.rds')