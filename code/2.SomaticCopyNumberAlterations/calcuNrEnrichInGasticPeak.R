library(dplyr)


nrGisticPeakQvalue <- readRDS(file = '/data/nrGisticPeakQvalue.rds')
ogGisticPeakQvalue <- readRDS(file = '/data/ogGisticPeakQvalue.rds')

############
nrPeakQvalue <- as.data.frame(nrGisticPeakQvalue)
nrPeakQvalue <- subset(nrPeakQvalue, !(Disease %in% c('COADREAD', 'GBMLGG', 'KIPAN', 'STES')))
nrPeakQvalue <- nrPeakQvalue %>% select(Approved.symbol, Disease, Direction) %>% mutate(Label = 'nrGene')


ogPeakQvalue <- as.data.frame(ogGisticPeakQvalue)
ogPeakQvalue <- subset(ogPeakQvalue, !(Disease %in% c('COADREAD', 'GBMLGG', 'KIPAN', 'STES')))

ogPeakQvalue <- ogPeakQvalue %>% select(gene_name, Disease, Direction) %>% mutate(Label = 'orGene') %>% 
  dplyr::rename(Approved.symbol = gene_name) 


gPeak <- rbind.data.frame(nrPeakQvalue, ogPeakQvalue)

############

geneStat <- data.frame(Label = c('nrGene', 'orGene'), gNum = c(111, 18638))
gPeak <- gPeak %>% group_by(Disease, Direction, Label) %>% count() %>% left_join(geneStat, by = 'Label')

gPeak <- gPeak %>% mutate(gNum = gNum - n)
fisherT <- function(data) {
  if(nrow(data) > 1){
    
    res <- fisher.test(data[, c('n', 'gNum')])
    
  }else{
    
    res <- fisher.test(rbind.data.frame(data.frame(n = 0, gNum = 111), data[, c('n', 'gNum')]))
    
  }
  
  return(data.frame(pValue = res$p.value, OR = res$estimate))
}



enrichResideInPeak <- gPeak %>% group_by(Disease, Direction) %>% do(fisherT(.))


saveRDS(enrichResideInPeak, file = '/data/enrichResideInPeak.rds')
