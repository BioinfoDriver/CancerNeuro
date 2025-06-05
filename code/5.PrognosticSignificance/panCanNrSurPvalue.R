
library(tibble)
library(survival)

####
tcgaPanCanCliData <- readRDS('/data/tcgaPanCanCliData.rds')
load(file = '/data/panCanGeneExpData.RData')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


neurotransmitterReceptors <- neurotransmitterReceptors %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID))
panCanNrExp <- panCanTurGeneExp[neurotransmitterReceptors$NCBI.Gene.ID, ]
rownames(panCanNrExp) <- neurotransmitterReceptors$Approved.symbol
colnames(panCanNrExp) <- substr(colnames(panCanNrExp), 1, 12)


comSams <- intersect(colnames(panCanNrExp), rownames(tcgaPanCanCliData))
panCanNrExp <- panCanNrExp[, comSams]
tcgaPanCanCliData <- tcgaPanCanCliData[comSams, ]

tcgaPanCanCliData <- tcgaPanCanCliData %>% rownames_to_column(var = 'patientId')
panCanNrExp <- as.data.frame(t(panCanNrExp)) %>% rownames_to_column(var = 'patientId')


panCanNrExpCliData <- merge.data.frame(panCanNrExp, tcgaPanCanCliData, by = 'patientId')

##############OS
osCanTypes <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
'LIHC','LUAD','LUSC','MESO','OV','PAAD','PRAD','READ','SARC','SKCM','STAD','THCA','UCEC','UCS','UVM')

nrGenes <- neurotransmitterReceptors$Approved.symbol

panCanSurPvalue <- lapply(osCanTypes, function(cancer){
  
  # canNrExpCliData <- subset(panCanNrExpCliData, cancer_type == 'ACC' & !is.na(os) & !is.na(os_time))
  
  canNrExpCliData <- subset(panCanNrExpCliData, cancer_type == cancer & !is.na(os) & !is.na(os_time))
  
  surPvalue <- lapply(nrGenes, function(gene){
    
    # gene <- 'HTR2B'
    
    naSum <- sum(is.na(canNrExpCliData[, gene]) | canNrExpCliData[, gene] == 0, na.rm = TRUE)
    
    if(naSum/nrow(canNrExpCliData) > 0.2){
      
      pValue <- c(NA, NA, NA)
      
    }else{
      
      
      nrExpCliData <- canNrExpCliData[, c('os_time', 'os', gene)]
      colnames(nrExpCliData) <- c('time', 'event', 'geneExp')
      
      nrExpCliData <- subset(nrExpCliData, !is.na(geneExp)) #  & geneExp > 0
      nrExpCliData <- nrExpCliData %>% mutate(expGroup = geneExp > median(geneExp))
      
      
      univ.model <- coxph(Surv(time, event)~geneExp, data = nrExpCliData)
      univ.p <- summary(univ.model)$coefficients[1, 5]
      hr <- summary(univ.model)$coefficients[1, 2]
      
      
      lr.test <- survdiff(Surv(time, event)~expGroup, data = nrExpCliData)
      lr.p <- 1 - pchisq(lr.test$chisq, length(lr.test$n) - 1)
      
      pValue <- c(hr, univ.p, lr.p)
      
    }
    
    names(pValue) <- c('HR', 'coxPvalue', 'logrankPvalue')
    
    return(pValue)
  })
  
  surPvalue <- as.data.frame(do.call(rbind, surPvalue))
  
  surPvalue$disease <- cancer
  surPvalue$geneSymbol <- nrGenes
  
  
  return(surPvalue)
})

panCanSurPvalue <- do.call(rbind.data.frame, panCanSurPvalue)


saveRDS(panCanSurPvalue, file = '/data/panCanNrSurPvalue.rds')


