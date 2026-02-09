

nReceptors <- readRDS(file = 'D:/CancerNeuroscience/Github/data/neurotransmitterReceptors.rds')
nReceptors <- nReceptors %>% mutate(NCBI.Gene.ID == as.character(NCBI.Gene.ID))

targetScan <- read.table(file = 'D:/CancerNeuroscience/data/TargetScan/Predicted_Targets_Context_Scores.default_predictions.txt', 
                         header = T, sep = '\t', stringsAsFactors = F)

nROfmiRNA <- subset(targetScan, Gene.Tax.ID == 9606) %>% mutate(Gene.ID = sub("\\..*", "", Gene.ID)) %>% 
  subset(Gene.ID %in% nReceptors$Ensembl.gene.ID) %>% select(Gene.Symbol, miRNA) %>% distinct()

load(file = 'D:/CancerNeuroscience/Github/data/expOfmiRNAs.RData')
# expOfmiRNAsStringent, expOfmiRNAsLoose
load(file = 'D:/CancerNeuroscience/Github/data/expOfNrGenes.RData')
# expOfNrGenesLoose, expOfNrGenesStringent


expOfNrGenesStringent <- lapply(expOfNrGenesStringent, function(expOfNrGenes){
  
  expOfNrGenes <- subset(nReceptors, NCBI.Gene.ID %in% expOfNrGenes)$Approved.symbol
  
  return(expOfNrGenes)
})


nROfmiRNAByDis <- sapply(names(expOfmiRNAsLoose), function(disease){
  
  miRNATarget <- subset(nROfmiRNA, Gene.Symbol %in% expOfNrGenesStringent[[disease]] & 
           miRNA %in% expOfmiRNAsStringent[[disease]])
    
  miRNATarget <- split(miRNATarget$miRNA, f = miRNATarget$Gene.Symbol)
  
  return(miRNATarget)
})


saveRDS(nROfmiRNAByDis, file = 'D:/CancerNeuroscience/Github/data/nROfmiRNAByDisease.rds')


