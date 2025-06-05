library('dplyr')
library('tibble')
library('cometExactTest')

###################
durgGene <- read.csv(file = '/data/oncokb_biomarker_drug_associations.tsv', header = T, sep = '\t', stringsAsFactors = F)
durgGene <- durgGene %>% subset(Level %in% c('1', '2', '3', '4'))
durgGeneList <- unique(durgGene$Gene)


neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
nrGeneList <- neurotransmitterReceptors$Approved.symbol

###########
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')
mc3MutData <- mc3MutData %>% mutate(score = 
                                      ifelse(Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',  'In_Frame_Ins', 
                                                                           'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'), 1, 0))

mc3MutData <- subset(mc3MutData, score == 1)
mc3MutData <- mc3MutData[!duplicated(mc3MutData[, c('Hugo_Symbol', 'SAMPLE_BARCODE')]), ]

mutData <- mc3MutData %>% reshape2::dcast(SAMPLE_BARCODE ~ Hugo_Symbol, value.var = 'score', fill = 0) %>% 
  column_to_rownames(var = "SAMPLE_BARCODE")


mutData <- mutData[, intersect(colnames(mutData), c(durgGeneList, nrGeneList))]
colnames(mutData) <- paste0(colnames(mutData), '_MUT')

###########
load(file = '/data/panCanNrCopyNumData.RData')
# nRCopyNumData, copyNumData

copyNumData <- copyNumData[intersect(rownames(copyNumData), c(durgGeneList, nrGeneList)), ]

delData <- apply(copyNumData, 2, function(x) ifelse(x == -2, 1, 0))
delData <- as.data.frame(t(delData))
colnames(delData) <- paste0(colnames(delData), '_DEL')


ampData <- apply(copyNumData, 2, function(x) ifelse(x == 2, 1, 0))
ampData <-  as.data.frame(t(ampData))
colnames(ampData) <- paste0(colnames(ampData), '_AMP')


###########
comSams <- Reduce(intersect, list(rownames(mutData), rownames(delData), rownames(ampData)))
mutData <- mutData[comSams, ]
delData <- delData[comSams, ]
ampData <- ampData[comSams, ]
altData <- cbind.data.frame(mutData, delData, ampData)

###################
nrGeneList <- c(paste0(nrGeneList, '_MUT'), paste0(nrGeneList, '_AMP'), paste0(nrGeneList, '_DEL'))
durgGeneList <- c(paste0(durgGeneList, '_MUT'), paste0(durgGeneList, '_AMP'), paste0(durgGeneList, '_DEL'))


tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
diseaseList <- unique(tcgaPanCanSamples$DISEASE)

# diseaseList <- 'LUSC'
panCanMeRes <- lapply(diseaseList, function(disease){
  
  # disease <- 'LUSC'
  diseaseSams <- subset(tcgaPanCanSamples, DISEASE == disease)$SAMPLE_BARCODE
  altMatrix <- altData[intersect(diseaseSams, rownames(altData)), ]
  
  altFre <- colSums(altMatrix)/nrow(altMatrix)
  altMatrix <- altMatrix[, names(altFre[altFre > 0.05])]
  
  nrNames <- intersect(nrGeneList, colnames(altMatrix))
  dgNames  <- intersect(durgGeneList, colnames(altMatrix))
  
  if(length(nrNames) > 1 & length(dgNames)){
    
    combs <- expand.grid(nrNames, dgNames)
    
    meRes <- lapply(seq(nrow(combs)), function(index){
      
      mutStat <- table(altMatrix[, combs[index, 1]], altMatrix[, combs[index, 2]])
      
      pValue <- comet_exact_test(c(mutStat[1, 1], mutStat[1, 2], mutStat[2, 1], mutStat[2, 2]))
      
      res <- data.frame(Disease = disease, nrGene = combs[index, 1], drGene = combs[index, 2], pValue = pValue)
      
    })
    
    meRes <- do.call(rbind, meRes)
    meRes$qValue <- p.adjust(meRes$pValue, method = 'fdr')
    
    return(meRes)
    
  }else{
    
    return(NULL)
    
  }
  
})

panCanMeRes <- do.call(rbind, panCanMeRes)

save(panCanMeRes, altData, file = '/data/panCanceMeRes.RDS')


