
library(dplyr)
library(tidyr)
library(tibble)

#################
nrReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')
durgTargetGenes <- read.csv(file = '/data/oncokb_biomarker_drug_associations.tsv', header = T, sep = '\t', stringsAsFactors = F)



#################
mc3MutData <- subset(mc3MutData, Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',  'In_Frame_Ins', 
                                                               'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))

mutData <- mc3MutData %>% distinct(SAMPLE_BARCODE, Hugo_Symbol, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = SAMPLE_BARCODE, names_from = Hugo_Symbol, values_from = Variant_Classification) %>% 
  column_to_rownames(var = "SAMPLE_BARCODE") %>% mutate_all(~ifelse(is.na(.), FALSE, TRUE)) 


rhoMut <- mutData[, intersect(colnames(mutData), nrReceptors$Approved.symbol)] %>% rename_with(~paste0(., "_MUT"), everything())

durgTargetMut <- mutData[, intersect(colnames(mutData), durgTargetGenes$Gene)] %>% rename_with(~paste0(., "_MUT"), everything())


#################
load('/data/panCanNrCopyNumData.RData')

rhoCnAmp <- as.data.frame(t(nRCopyNumData)) %>% mutate_all(~ifelse(. == 2, TRUE, FALSE)) %>% rename_with(~paste0(., "_AMP"), everything())
rhoCnDel <- as.data.frame(t(nRCopyNumData)) %>% mutate_all(~ifelse(. == -2, TRUE, FALSE)) %>% rename_with(~paste0(., "_DEL"), everything())

indexes <- colSums(rhoCnAmp) >= colSums(rhoCnDel)
rhoCnAmp <- rhoCnAmp[, indexes]
rhoCnDel <- rhoCnDel[, !indexes]



durgTargetCnAmp <- as.data.frame(t(copyNumData[intersect(rownames(copyNumData),durgTargetGenes$Gene ), ])) %>% 
  mutate_all(~ifelse(. == 2, TRUE, FALSE)) %>% rename_with(~paste0(., "_AMP"), everything())
durgTargetCnDel <- as.data.frame(t(copyNumData[intersect(rownames(copyNumData),durgTargetGenes$Gene ), ])) %>% 
  mutate_all(~ifelse(. == -2, TRUE, FALSE)) %>% rename_with(~paste0(., "_DEL"), everything())

indexes <- colSums(durgTargetCnAmp) >= colSums(durgTargetCnDel)
durgTargetCnAmp <- durgTargetCnAmp[, indexes]
durgTargetCnDel <- durgTargetCnDel[, !indexes]


#################
gam <- Reduce(function(x, y){
  res <- merge(x, y, by = "row.names", all = FALSE) %>% column_to_rownames(var = "Row.names")
  return(res)
}, list(rhoMut, durgTargetMut, rhoCnAmp, rhoCnDel, durgTargetCnAmp, durgTargetCnDel))


alteration.class <- setNames(gsub(".*_", "", colnames(gam)), colnames(gam))


########
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

tcgaPanCanSamples <- tcgaPanCanSamples %>% mutate(sample.class = paste(DISEASE, SUBTYPE, sep = '_')) %>% 
  mutate(sample.class = gsub('_Not_Applicable', '', sample.class)) %>% subset(SAMPLE_BARCODE %in% rownames(gam))

########
sample.class <- setNames(tcgaPanCanSamples$sample.class, tcgaPanCanSamples$SAMPLE_BARCODE)

save(gam, alteration.class, sample.class, file = '/data/dataForSelectCal.RData')

