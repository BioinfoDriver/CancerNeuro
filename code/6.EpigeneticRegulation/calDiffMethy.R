library('dplyr')
library('tibble')

##############
load(file = '/data/panCanMethyData.RData')
# panCanTurMethy, panCanPairdTurMethy, panCanPairdNormMethy

anno450k <- readRDS(file = '/data/methAnno450k.rds')
nReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')

nrAnno450k <- subset(anno450k, RefGene_Name %in% nReceptors$Approved.symbol)
nrpromoterAnno <- subset(nrAnno450k, RefGene_Group %in% c("5'UTR", "TSS1500", "TSS200") | Enhancer == 'TRUE')
nrpromoterAnno <- nrpromoterAnno[, c('chr', 'pos', 'strand', 'Name', 'RefGene_Name', 'Enhancer', 'RefGene_Group')] %>% distinct()


panCanPairdTurMethy <- rownames_to_column(panCanPairdTurMethy, var = 'Name')
panCanPairdNormMethy <- rownames_to_column(panCanPairdNormMethy, var = 'Name')


panCanPairdTurMethy <- panCanPairdTurMethy %>% inner_join(nrpromoterAnno, by = 'Name')
panCanPairdNormMethy <- panCanPairdNormMethy %>% inner_join(nrpromoterAnno, by = 'Name')


###############
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


pairdSams <- substr(colnames(panCanPairdTurMethy), 1, 12)
pairdSams <- tcgaPanCanSamples %>% subset(PATIENT_BARCODE %in% pairdSams)


samSta <- pairdSams %>% group_by(DISEASE) %>% count() %>% subset(n >=10)
pairdSams <- subset(pairdSams, DISEASE %in% samSta$DISEASE)


tumorSams <- pairdSams %>% mutate(samBarcode = SAMPLE_BARCODE, samType = 'Tumor')
normalSams <- pairdSams %>% mutate(samBarcode = paste0(PATIENT_BARCODE, '-11'), samType = 'Normal')


pairdSams <- rbind.data.frame(tumorSams, normalSams)
pairdSams <- split.data.frame(pairdSams, f = pairdSams$DISEASE)


# "BLCA" "BRCA" "COAD" "ESCA" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC"

# differently methylation
diseases <- names(pairdSams)

diffMethy <- lapply(diseases, function(disease){
  
  diseaseSams <- pairdSams[[disease]]
  
  tM <- panCanPairdTurMethy[, subset(diseaseSams, samType == 'Tumor')$samBarcode]
  nM <- panCanPairdNormMethy[, subset(diseaseSams, samType == 'Normal')$samBarcode]

  
  diffM <- lapply(seq(nrow(tM)), function(index){
    
    tnv <- data.frame(tMv = as.numeric(tM[index, ]), nMv = as.numeric(nM[index, ]))
    ftnv <- subset(tnv, !is.na(tMv) & !is.na(nMv))
    
    if(nrow(ftnv)/nrow(tnv) >= 0.9){
      
      logFC <- log2(mean(ftnv$tMv)/mean(ftnv$nMv))
      mDiff <- mean(ftnv$tMv) - mean(ftnv$nMv)
      pValue <- wilcox.test(Pair(tMv, nMv) ~ 1, data = ftnv)$p.value
      
    }else{
      
      logFC <- NA
      mDiff <- NA
      pValue <- NA
      
    }
    
    return(data.frame(logFC, mDiff, pValue))
  })

  diffM <- do.call(rbind, diffM)
  diffM <- cbind.data.frame(panCanPairdTurMethy[, c('chr', 'pos', 'strand', 'Name', 'RefGene_Name', 'Enhancer', 'RefGene_Group')], diffM)
  
  diffM <- diffM %>% subset(!is.na(logFC)) %>% mutate(adjPVal = p.adjust(pValue, method = 'fdr'), Disease = disease)
  return(diffM)
})

diffMethy <- do.call(rbind.data.frame, diffMethy)


# Determine the threshold
diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0.585 & adjPVal < 0.05, 'Hyper', ifelse(logFC < -0.585 & adjPVal < 0.05, 'Hypo', 'Neutral')))

mDiffDensityPlot <- ggplot(subset(diffMethy, Mstatus != 'Neutral'), aes(x = mDiff)) +geom_density(fill = "skyblue") + 
  labs(x = "Average methylation differences", y = "Density") + scale_x_continuous(breaks = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)) + 
  annotate("text", x = -0.05, y = 0.5, label = "-0.1", color = 'blue') + 
  annotate("text", x = 0.15, y = 3.8, label = "0.1", color = 'blue') + 
  geom_vline(xintercept = -0.1, linetype = 'dashed', color = 'red', size = 0.5) + 
  geom_vline(xintercept = 0.1, linetype = 'dashed', color = 'red', size = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))

ggsave(filename = '/result/section4/mDiffDensityPlot.pdf', 
       plot = mDiffDensityPlot, width = 7, height = 7, units = 'cm')


diffMethy <- diffMethy %>% mutate(Mstatus = ifelse(logFC > 0.585 & adjPVal < 0.05 & mDiff > 0.1, 'Hyper', 
                                                   ifelse(logFC < -0.585 & adjPVal < 0.05 & mDiff < -0.1, 'Hypo', 'Neutral')))
colnames(diffMethy)[c(5, 12)] <- c('Symbol', 'DISEASE')


# Filter out inconsistent entries
diffMethy <- lapply(split.data.frame(diffMethy, f = ~ DISEASE + Symbol), function(methyStatus){
  
  if(any(methyStatus$Mstatus == 'Hyper') & any(methyStatus$Mstatus == 'Hypo')){
    
    return(NULL)
    
  }else{
    
    return(methyStatus)
  }

})

diffMethy <- Filter(Negate(is.null), diffMethy)
diffMethy <- do.call(rbind.data.frame, diffMethy)
diffMethy <- diffMethy %>% remove_rownames()

saveRDS(diffMethy, file = '/data/panCanDiffMethy.rds')
