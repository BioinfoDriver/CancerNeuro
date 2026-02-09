
nROfmiRNAByDis <- readRDS(file = '/data/nROfmiRNAByDisease.rds')

nROfmiRNAByDis <- sapply(nROfmiRNAByDis, function(nROfmiRNA){
  
  nrRegulators <- sapply(nROfmiRNA, function(nrRegulator){
    
    nrRegulator <- gsub(pattern = '-', replacement = '\\.', nrRegulator)
    
    return(nrRegulator)
  })
  
  return(nrRegulators)
})


nReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
nReceptors <- nReceptors %>% mutate(NCBI.Gene.ID == as.character(NCBI.Gene.ID))


load('/data/panCanMiRnaExpData.RData')
rownames(panCanTurMiRnaExp) <- gsub(pattern = '-', replacement = '\\.', rownames(panCanTurMiRnaExp))

load('/data/panCanGeneExpData.RData')
panCanTurNrExp <- panCanTurGeneExp[nReceptors$NCBI.Gene.ID, ]
rownames(panCanTurNrExp) <- nReceptors$Approved.symbol


tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
tcgaPanCanSamples <- split.data.frame(tcgaPanCanSamples, f = ~DISEASE)


# Fitting Linear Models
nROfmiRNAflm <- lapply(names(nROfmiRNAByDis), function(disease){
  
  # disease <- 'ACC'
  nrExp <- panCanTurNrExp %>% select(any_of(tcgaPanCanSamples[[disease]]$SAMPLE_BARCODE))%>%
    filter(rownames(.) %in% names(nROfmiRNAByDis[[disease]]))
  miRNAExp <- panCanTurMiRnaExp %>% select(any_of(tcgaPanCanSamples[[disease]]$SAMPLE_BARCODE))%>%
    filter(rownames(.) %in% unique(unlist(nROfmiRNAByDis[[disease]])))
  
  commonSams <- intersect(names(nrExp), names(miRNAExp))
  nRmiRNAExpData <- rbind(nrExp[, commonSams], miRNAExp[, commonSams])
  
  nRmiRNAExpData[is.na(nRmiRNAExpData)] <- 0
  nRmiRNAExpData <- as.data.frame(t(nRmiRNAExpData))

  
  coefPvalues <- lapply(names(nROfmiRNAByDis[[disease]]), function(nrGene){
    
    # nrGene <- "ADORA1"
    f = as.formula(paste0(nrGene, ' ~ ', paste0(nROfmiRNAByDis[[disease]][[nrGene]], collapse = '+')))
    model <- lm(formula = f, data = nRmiRNAExpData)
    model <- summary(model)
    
    coefPvalue <- as.data.frame(model$coefficients)[nROfmiRNAByDis[[disease]][[nrGene]], ]
    
    coefPvalue <- coefPvalue %>% rownames_to_column(var = 'miRNA') %>% mutate(DISEASE = disease, Gene = nrGene)
    return(coefPvalue)
  })
  
  coefPvalues <- do.call(rbind, coefPvalues)
  
  return(coefPvalues)
})


nROfmiRNAflm <- do.call(rbind, nROfmiRNAflm)

colnames(nROfmiRNAflm) <- c('miRNA', 'Estimate', 'StdError', 'tValue', 'pValue', 'Disease', 'nRgene')
nROfmiRNAflm <- nROfmiRNAflm %>% mutate(miRNA = gsub(pattern = '\\.', replacement = '-', miRNA))

saveRDS(nROfmiRNAflm, file = '/data/nROfmiRNARegulatorFLM.rds')



