library('dplyr')
library('ggplot2')
library('readxl')
library('GSVA')
library('ggpubr')
############
setwd('/data')
load(file = 'panCanceMeRes.RDS')
# panCanMeRes, altData, 

tcgaPanCanSamples <- readRDS(file = 'tcgaPanCanSamples.rds')


load(file = 'panCanGeneExpData.RData')
# geneInfo, panCanTurGeneExp, panCanPairdTurGeneExp, panCanPairdNormGeneExp


# Activation of cAMP-Dependent PKA
cAMP_Pathway <- read_excel(path = 'cAMP_Pathway.xlsx', col_names = FALSE)
cAMP_Pathway <- c(cAMP_Pathway[, 1, T], cAMP_Pathway[, 2, T], cAMP_Pathway[, 3, T], 
                  cAMP_Pathway[, 4, T], cAMP_Pathway[, 5, T], cAMP_Pathway[, 6, T], cAMP_Pathway[, 7, T], cAMP_Pathway[, 8, T])


cAMP_Pathway <- geneInfo$geneID[na.omit(match(cAMP_Pathway, geneInfo$geneName))]

# proliferating score——median score
cAMPMedianScore <- apply(panCanTurGeneExp[cAMP_Pathway, ], 2, median, na.rm = T)	
cAMPScore <- cbind(altData, prolif = cAMPMedianScore[rownames(altData)])
cAMPScore$DISEASE <- tcgaPanCanSamples$DISEASE[match(rownames(cAMPScore), tcgaPanCanSamples$SAMPLE_BARCODE)]



# Tumor cAMP score compare
wilcox.test(prolif~ADRB3_AMP, cAMPScore, alternative = 'two.sided')$p.value # 3.402595e-23
wilcox.test(prolif~ADRB3_AMP, subset(cAMPScore, DISEASE == 'BRCA'), alternative = 'two.sided')$p.value
wilcox.test(prolif~ADRB3_AMP, subset(cAMPScore, DISEASE == 'ESCA'), alternative = 'two.sided')$p.value
wilcox.test(prolif~ADRB3_AMP, subset(cAMPScore, DISEASE == 'HNSC'), alternative = 'two.sided')$p.value



HNSCcAMPScore <- subset(cAMPScore, DISEASE == 'HNSC') %>% mutate(label = paste(ADRB3_AMP, ATR_AMP, sep = '_'))
ESCAcAMPScore <- subset(cAMPScore, DISEASE == 'ESCA') %>% mutate(label = paste(ADRB3_AMP, CDK12_AMP, sep = '_'))
BRCAcAMPScore <- subset(cAMPScore, DISEASE == 'BRCA') %>% mutate(label = paste(ADRB3_AMP, TP53_MUT, sep = '_'))


ggplot(data = cAMPScore, aes(as.factor(ADRB3_AMP), prolif)) +
  geom_boxplot(aes(fill = ADRB3_AMP), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  xlab('Subtype') + ylab('Proliferating index') + guides(fill = guide_legend(title= 'Subtype')) + 
  theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')





# Tumor proliferating gene
pcnaSig <- read.csv(file = 'PCNA_Signature.txt', header = TRUE,sep = '\t',stringsAsFactors =FALSE)
pcnaSig <- intersect(pcnaSig$Entrez.ID, rownames(panCanTurGeneExp))


# proliferating score——median score
pcnaMedianScore <- apply(panCanTurGeneExp[pcnaSig, ], 2, median, na.rm = T)	


#  proliferating score——GSVA score
pcnaGsvaScore <- gsva(expr=as.matrix(panCanTurGeneExp), gset.idx.list=list(pcna=pcnaSig), method="gsva")
pcnaGsvaScore <- as.data.frame(t(pcnaGsvaScore))

pcnaSsgseaScore <- gsva(expr=as.matrix(panCanTurGeneExp), gset.idx.list=list(pcna=pcnaSig), method="ssgsea")
pcnaSsgseaScore <- as.data.frame(t(pcnaSsgseaScore))



# prolifScore <- cbind(altData, prolif = pcnaMedianScore[rownames(altData)])
# prolifScore <- cbind(altData, prolif = pcnaGsvaScore[rownames(altData), ])
prolifScore <- cbind(altData, prolif = pcnaSsgseaScore[rownames(altData), ])
prolifScore$DISEASE <- tcgaPanCanSamples$DISEASE[match(rownames(prolifScore), tcgaPanCanSamples$SAMPLE_BARCODE)]


# Tumor proliferating score compare
wilcox.test(prolif~ADRB3_AMP, prolifScore, alternative = 'two.sided')$p.value # 3.402595e-23
wilcox.test(prolif~ADRB3_AMP, subset(prolifScore, DISEASE == 'BRCA'), alternative = 'two.sided')$p.value
wilcox.test(prolif~ADRB3_AMP, subset(prolifScore, DISEASE == 'ESCA'), alternative = 'two.sided')$p.value
wilcox.test(prolif~ADRB3_AMP, subset(prolifScore, DISEASE == 'HNSC'), alternative = 'two.sided')$p.value


HNSCprolifScore <- subset(prolifScore, DISEASE == 'HNSC') %>% mutate(label = paste(ADRB3_AMP, ATR_AMP, sep = '_'))
ESCAprolifScore <- subset(prolifScore, DISEASE == 'ESCA') %>% mutate(label = paste(ADRB3_AMP, CDK12_AMP, sep = '_'))
BRCAprolifScore <- subset(prolifScore, DISEASE == 'BRCA') %>% mutate(label = paste(ADRB3_AMP, TP53_MUT, sep = '_'))


ggplot(data = HNSCprolifScore, aes(label, prolif)) +
  geom_boxplot(aes(fill = label), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  xlab('Subtype') + ylab('Proliferating index') + guides(fill = guide_legend(title= 'Subtype')) + 
  theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')


