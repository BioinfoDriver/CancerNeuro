
library('dplyr')
library('tibble')


load(file = '/data/randomGenesList.RData')
load(file = '/data/neuroRecepCnStat.RData')
load(file = '/data/getPanCanNrCopyNumData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


# by disease
rgCnStatByDis <- lapply(randomGenesList, function(randomGenes){
  
  # randomGenes <- randomGenesList[[1]]
  rgCopyNumData <- copyNumData[intersect(rownames(copyNumData), randomGenes$gene_name), ]
  
  
  rgCopyNumData <- rownames_to_column(rgCopyNumData, var = 'Approved.symbol') %>% 
    reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')
  
  rgCopyNumData <- rgCopyNumData %>% left_join(tcgaPanCanSamples, by = join_by(SAMPLE_BARCODE))
  
  
  samsStat <- rgCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, cnStaus) %>% summarise(n = n()) %>% 
    group_by(DISEASE)  %>% summarize(numberOfPatient = n_distinct(SAMPLE_BARCODE))
  
  ##########
  rgCnStatByDis <- rgCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, cnStaus) %>% count() %>% 
    group_by(DISEASE, cnStaus) %>% summarize(numberOfPatientAlt = n())
  
  rgCnStatByDis <- rgCnStatByDis %>% subset(cnStaus != 0) %>% left_join(samsStat, by = join_by(DISEASE)) %>% 
    mutate(altFre = numberOfPatientAlt/numberOfPatient) %>% select(DISEASE, cnStaus, altFre)
  
  return(rgCnStatByDis)
})

rgCnStatByDis <- Reduce(function(x, y) merge(x = x, y = y, by = c('DISEASE', 'cnStaus'), all = TRUE), rgCnStatByDis)


colnames(rgCnStatByDis) <- c('DISEASE', 'cnStaus', paste0('rgRes', 1:1000))

comResByDis <- merge(nRCnStatByDis, rgCnStatByDis, by = c('DISEASE', 'cnStaus'), all = TRUE)
comResByDis <- comResByDis %>% arrange(desc(DISEASE), desc(cnStaus))
comResByDis <- split.data.frame(comResByDis, f =  ~ DISEASE + cnStaus)


comResByDisPvalue <- lapply(comResByDis, function(comRes){
  
  comRes$pValue <- 1- sum(comRes[1, 6] >= comRes[1, 6:1005], na.rm = TRUE)/1000

  return(comRes)
})
comResByDisPvalue <- do.call(rbind, comResByDisPvalue)



# by disease/group
rgCnStatListByDisGroup <- lapply(randomGenesListByGroup, function(rgGenesByGroup){
  
  rgCnStatByDisGroup <- lapply(rgGenesByGroup, function(rgGenes){
    
    rgCopyNumData <- copyNumData[intersect(rownames(copyNumData), rgGenes$gene_name), ]
    
    
    rgCopyNumData <- rownames_to_column(rgCopyNumData, var = 'Approved.symbol') %>% 
      reshape2::melt(id.vars = 'Approved.symbol', variable.name = 'SAMPLE_BARCODE', value.name = 'cnStaus')
    
    rgCopyNumData <- rgCopyNumData %>% left_join(tcgaPanCanSamples, by = join_by(SAMPLE_BARCODE))
    
    
    samsStat <- rgCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, cnStaus) %>% summarise(n = n()) %>% 
      group_by(DISEASE)  %>% summarize(numberOfPatient = n_distinct(SAMPLE_BARCODE))
    
    ##########
    rgCnStatByDis <- rgCopyNumData %>% group_by(DISEASE, SAMPLE_BARCODE, cnStaus) %>% count() %>% 
      group_by(DISEASE, cnStaus) %>% summarize(numberOfPatientAlt = n())
    
    rgCnStatByDis <- rgCnStatByDis %>% subset(cnStaus != 0) %>% left_join(samsStat, by = join_by(DISEASE)) %>% 
      mutate(altFre = numberOfPatientAlt/numberOfPatient) %>% select(DISEASE, cnStaus, altFre)
    
    return(rgCnStatByDis)
    
  })
  
  rgCnStatByDisGroup <- Reduce(function(x, y) merge(x = x, y = y, by = c('DISEASE', 'cnStaus'), all = TRUE), rgCnStatByDisGroup)
  
  colnames(rgCnStatByDisGroup) <- c('DISEASE', 'cnStaus', paste0('rgRes', 1:1000))
  
  return(rgCnStatByDisGroup)
})


groupNames <- names(table(nRCnStatByGroupDis$classesOfNeuroreceptors))

comResByDisGroupPvalue <- lapply(seq(length(groupNames)), function(i){
  
  gName <- groupNames[i]
  nrCnStatByDisGroup <- subset(nRCnStatByGroupDis, classesOfNeuroreceptors == gName)
  
  
  comResByDisGroup <- merge(nrCnStatByDisGroup, rgCnStatListByDisGroup[[i]], by = c('DISEASE', 'cnStaus'), all = TRUE)
  comResByDisGroup[is.na(comResByDisGroup)] <- 0
  
  comResByDisGroupPvalue <- sapply(seq(nrow(comResByDisGroup)), function(i){
    
    res <- sum(comResByDisGroup[i, 6] >= comResByDisGroup[i, 6:1005], na.rm = TRUE)
    return(res)
    
  })
  names(comResByDisGroupPvalue) <- paste(comResByDisGroup$DISEASE, comResByDisGroup$cnStaus, sep = '_')  
  comResByDisGroupPvalue <- 1 - comResByDisGroupPvalue/1000
  
  return(comResByDisGroupPvalue)
})

comResByDisGroupPvalue <- do.call(cbind, comResByDisGroupPvalue)
colnames(comResByDisGroupPvalue) <- groupNames



save(comResByDisPvalue, comResByDisGroupPvalue, file = '/data/neuroRecepCnComToRandom.RData')


