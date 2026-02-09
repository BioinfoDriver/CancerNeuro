library('dplyr')
library('tibble')


load('/data/panCanMiRnaExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

tcgaPanCanSamples <- split.data.frame(tcgaPanCanSamples, f = ~DISEASE)


expOfmiRNAsStringent <- sapply(names(tcgaPanCanSamples), function(disease){
  
  expData <- panCanTurMiRnaExp %>% select(any_of(tcgaPanCanSamples[[disease]]$SAMPLE_BARCODE))
  
  filterIndex <- apply(expData, 1, function(x) sum(x == 0, na.rm = TRUE) + sum(is.na(x)))/ncol(expData) < 0.1
  
  expOfNrs <- rownames(expData)[filterIndex]
  
  return(expOfNrs)
})


expOfmiRNAsLoose <- sapply(names(tcgaPanCanSamples), function(disease){
  
  expData <- panCanTurMiRnaExp %>% select(any_of(tcgaPanCanSamples[[disease]]$SAMPLE_BARCODE))
  
  filterIndex <- apply(expData, 1, function(x) sum(x == 0, na.rm = TRUE) + sum(is.na(x)))/ncol(expData) < 0.2
  
  expOfNrs <- rownames(expData)[filterIndex]
  
  return(expOfNrs)
})


save(expOfmiRNAsStringent, expOfmiRNAsLoose, file = '/data/expOfmiRNAs.RData')



