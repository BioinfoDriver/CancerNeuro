
library('limma')
library('dplyr')
library('tibble')

###############
load(file = 'D:/CancerNeuroscience/Github/data/panCanMiRnaExpData.RData')
tcgaPanCanSamples <- readRDS(file = 'D:/CancerNeuroscience/Github/data/tcgaPanCanSamples.rds')

load(file = 'D:/CancerNeuroscience/Github/data/expOfmiRNAs.RData')
# expOfmiRNAsStringent, expOfmiRNAsLoose

pairdSams <- substr(colnames(panCanPairdTurMiRnaExp), 1, 12)
pairdSams <- tcgaPanCanSamples %>% subset(PATIENT_BARCODE %in% pairdSams)

samSta <- pairdSams %>% group_by(DISEASE) %>% count() %>% subset(n >=10)
pairdSams <- subset(pairdSams, DISEASE %in% samSta$DISEASE)

tumorSams <- pairdSams %>% mutate(samBarcode = SAMPLE_BARCODE, samType = 'Tumor')
normalSams <- pairdSams %>% mutate(samBarcode = paste0(PATIENT_BARCODE, '-11'), samType = 'Normal')

pairdSams <- rbind.data.frame(tumorSams, normalSams)
pairdSams <- split.data.frame(pairdSams, f = pairdSams$DISEASE)


# sapply(pairdSams, nrow)/2
# BLCA BRCA ESCA HNSC KICH KIRC KIRP LIHC LUAD LUSC PRAD STAD THCA UCEC 
#   19   85   12   42   24   65   34   45   44   41   49   33   56   15 
# sum(sapply(pairdSams, nrow)/2)
# [1] 564


# differently expressed miRNA
diseases <- names(pairdSams)
mirDiffExp <- lapply(diseases, function(disease){
  
  expmiRNA <- expOfmiRNAsStringent[[disease]]
  
  diseaseSams <- pairdSams[[disease]]
  diseaseSams$PATIENT_BARCODE <- factor(diseaseSams$PATIENT_BARCODE)
  diseaseSams$samType <- factor(diseaseSams$samType, levels = c('Normal', 'Tumor'))
  
  design <- model.matrix(~PATIENT_BARCODE+samType, data = diseaseSams)
  
  
  eset <- cbind.data.frame(panCanPairdTurMiRnaExp[expmiRNA, subset(diseaseSams, samType == 'Tumor')$samBarcode], 
                           panCanPairdNormMiRnaExp[expmiRNA, subset(diseaseSams, samType == 'Normal')$samBarcode])
  
  # filter
  filterIndex <- apply(eset, 1, function(x) sum(x == 0, na.rm = TRUE) + sum(is.na(x)))/ncol(eset) < 0.1
  eset <- eset[filterIndex, ]
  
  
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  difRes <- topTable(fit, coef="samTypeTumor", adjust.method = "BH", n = Inf)
  
  difRes <- difRes %>% mutate(DISEASE = disease) %>% rownames_to_column(var = "geneID")
  return(difRes)
})

mirDiffExp <- do.call(rbind.data.frame, mirDiffExp)

saveRDS(mirDiffExp, file = 'D:/CancerNeuroscience/Github/data/panCanMiRnaDiffExp.rds')


