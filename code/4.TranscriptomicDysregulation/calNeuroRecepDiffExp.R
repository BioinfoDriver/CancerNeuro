
library('limma')
library('dplyr')
library('tibble')


############
setwd('/data')

load(file = 'panCanGeneExpData.RData')
tcgaPanCanSamples <- readRDS(file = 'tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = 'neurotransmitterReceptors.rds')


pairdSams <- substr(colnames(panCanPairdTurGeneExp), 1, 12)
pairdSams <- tcgaPanCanSamples %>% subset(PATIENT_BARCODE %in% pairdSams)


samSta <- pairdSams %>% group_by(DISEASE) %>% count() %>% subset(n >=10)
pairdSams <- subset(pairdSams, DISEASE %in% samSta$DISEASE)


tumorSams <- pairdSams %>% mutate(samBarcode = SAMPLE_BARCODE, samType = 'Tumor')
normalSams <- pairdSams %>% mutate(samBarcode = paste0(PATIENT_BARCODE, '-11'), samType = 'Normal')


pairdSams <- rbind.data.frame(tumorSams, normalSams)
pairdSams <- split.data.frame(pairdSams, f = pairdSams$DISEASE)


# differently expressed genes
neurotransmitterReceptors <- neurotransmitterReceptors %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID))
diseases <- names(pairdSams)

nRDiffExp <- lapply(diseases, function(disease){
  
  
  diseaseSams <- pairdSams[[disease]]
  
  
  diseaseSams$PATIENT_BARCODE <- factor(diseaseSams$PATIENT_BARCODE)
  diseaseSams$samType <- factor(diseaseSams$samType, levels = c('Normal', 'Tumor'))
                                         
  design <- model.matrix(~PATIENT_BARCODE+samType, data = diseaseSams)
  
  
  eset <- cbind.data.frame(panCanPairdTurGeneExp[, subset(diseaseSams, samType == 'Tumor')$samBarcode], 
                           panCanPairdNormGeneExp[, subset(diseaseSams, samType == 'Normal')$samBarcode])
  
  eset <- eset[neurotransmitterReceptors$NCBI.Gene.ID, ]
  
  # filter
  filterIndex <- apply(eset, 1, function(x) sum(x == 0, na.rm = TRUE) + sum(is.na(x)))/ncol(eset) > 0.2
  eset <- eset[!filterIndex, ]
  
  
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  difRes <- topTable(fit, coef="samTypeTumor", adjust.method = "BH", n = Inf)
  
  difRes <- difRes %>% mutate(DISEASE = disease) %>% rownames_to_column(var = "geneID")
  return(difRes)
})


nRDiffExp <- do.call(rbind.data.frame, nRDiffExp)
nRDiffExp <- merge(nRDiffExp, neurotransmitterReceptors[, c('NCBI.Gene.ID', 'Approved.symbol')], 
                   by.x = 'geneID', by.y = 'NCBI.Gene.ID')


saveRDS(nRDiffExp, file = 'panCanNrDiffExp.rds')


