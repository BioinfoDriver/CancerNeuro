
library('dplyr')

###################
durgGene <- read.csv(file = '/data/oncokb_biomarker_drug_associations.tsv', header = T, sep = '\t', stringsAsFactors = F)
durgGene <- durgGene %>% subset(Level %in% c('1', '2', '3', '4'))

load(file = '/data/panCanceMeRes.RDS')
# panCanMeRes, altData, 

panCanSigMeRes <- subset(panCanMeRes, qValue < 0.05)

panCanSigMeRes <- panCanSigMeRes %>% rename(nrGeneAltType = nrGene, drGeneAltType = drGene) %>% 
  mutate(nrGeneAltType = as.character(nrGeneAltType), drGeneAltType = as.character(drGeneAltType))

##########
nrGeneAltType <- do.call(rbind, strsplit(panCanSigMeRes$nrGeneAltType, split = "_"))
colnames(nrGeneAltType) <- c('nrNames', 'nrAltType')
drGeneAltType <- do.call(rbind, strsplit(panCanSigMeRes$drGeneAltType, split = "_"))  
colnames(drGeneAltType) <- c('drNames', 'drAltType')

##########
panCanSigMeRes <- cbind(panCanSigMeRes, nrGeneAltType, drGeneAltType)

panCanSigMeRes <- panCanSigMeRes %>% left_join(durgGene, by = join_by(drNames == Gene))


# write.table(panCanSigMeRes, file = '/result/section2/panCanSigMeRes.txt', row.names = T, sep = '\t', quote = F)


##########
mePaireds <- panCanSigMeRes %>% distinct(Disease, nrGeneAltType, drGeneAltType)

tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


sapply(seq(nrow(mePaireds)), function(index){
  
  mePaired <- mePaireds[index, ]
  
  sams <- subset(tcgaPanCanSamples, DISEASE == mePaired$Disease)$SAMPLE_BARCODE
  
  meMat <- altData[intersect(rownames(altData), sams), c(mePaired$nrGeneAltType, mePaired$drGeneAltType)]
  
  colnames(meMat) <- c('Gene1', 'Gene2')
  meMat <- meMat %>% arrange(desc(Gene1), desc(Gene2))
  
  altFre <- round(colSums(meMat)/nrow(meMat), 2)
  
  totalAltNum <- sum(meMat$Gene1 | meMat$Gene2)
  totalAltFre <- round(totalAltNum/nrow(meMat), 2)
  
  
  colnames(meMat) <- paste0(c(mePaired$nrGeneAltType, mePaired$drGeneAltType), '(', altFre, ')')
  
  
  
  title <- paste0('Coverage: ', totalAltFre, '(', totalAltNum, '/', nrow(meMat), ')', ' of samples of ', mePaired$Disease)
  
  
  setwd('/results/section2/mePlots/')
  pdf(file = paste0(mePaired$Disease, '_', mePaired$nrGeneAltType, '_', mePaired$drGeneAltType,'.pdf'))
  
   pheatmap::pheatmap(t(meMat), color = c('white', 'red'), border_color = 'gray', 
                      cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, main = title)
  
  dev.off()
  
  return(NULL)
})





