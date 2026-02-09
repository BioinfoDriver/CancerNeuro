
library('dplyr')
library('ggplot2')
library('ggExtra')
library('Hmisc')
library('RColorBrewer')
##########
tcgaPanCanSamples <- readRDS(file = 'D:/CancerNeuroscience/Github/data/tcgaPanCanSamples.rds')


load(file = 'D:/CancerNeuroscience/Github/data/panCanMiRnaExpData.RData')
# panCanTurMiRnaExp, panCanPairdTurMiRnaExp, panCanPairdNormMiRnaExp

load(file = 'D:/CancerNeuroscience/Github/data/panCanGeneExpData.RData')
# geneInfo, panCanTurGeneExp, panCanPairdTurGeneExp, panCanPairdNormGeneExp

nReceptors <- readRDS(file = 'D:/CancerNeuroscience/Github/data/neurotransmitterReceptors.rds')
nReceptors <- nReceptors %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID))

comSmas <- intersect(colnames(panCanTurMiRnaExp), colnames(panCanTurGeneExp))


turMiRnaExp <- panCanTurMiRnaExp[, comSmas]
turGeneExp <- panCanTurGeneExp[nReceptors$NCBI.Gene.ID, comSmas]
rownames(turGeneExp) <- nReceptors$Approved.symbol


tcgaPanCanSamples <- subset(tcgaPanCanSamples, SAMPLE_BARCODE %in% comSmas)


############
KIRCsamples <- subset(tcgaPanCanSamples, DISEASE == 'PRAD')
KIRCMiRnaExp <- turMiRnaExp['hsa-miR-339-5p', KIRCsamples$SAMPLE_BARCODE]
KIRCGeneExp <- turGeneExp['P2RY2', KIRCsamples$SAMPLE_BARCODE]


mirGeneExp <- as.data.frame(t(rbind.data.frame(KIRCMiRnaExp, KIRCGeneExp)))
colnames(mirGeneExp) <- c('miRNA', 'Target')


miRNATargetExpPlot <- ggplot(mirGeneExp, aes(x = miRNA, y = Target)) + geom_point() + geom_smooth(method = "loess", se = F) + 
  labs(x = "Expression of miR-339-5p", y = "Expression of P2RY2") + 
  theme_bw() + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))


############
comSmas <- Reduce(intersect, list(substr(colnames(panCanPairdTurMiRnaExp), 1, 12), substr(colnames(panCanPairdNormGeneExp), 1, 12),
substr(colnames(panCanPairdTurGeneExp), 1, 12), substr(colnames(panCanPairdNormMiRnaExp), 1, 12)))

sams <- subset(tcgaPanCanSamples, DISEASE == 'PRAD' & PATIENT_BARCODE %in% comSmas)



tMirExp <- panCanPairdTurMiRnaExp["hsa-miR-339-5p", paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tMirExp) <- substr(colnames(tMirExp), 1, 12)

tGeneExp <- panCanPairdTurGeneExp["5029", paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tGeneExp) <- substr(colnames(tGeneExp), 1, 12)

nMirExp <- panCanPairdNormMiRnaExp["hsa-miR-339-5p", paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nMirExp) <- substr(colnames(nMirExp), 1, 12)

nGeneExp <- panCanPairdNormGeneExp["5029", paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nGeneExp) <- substr(colnames(nGeneExp), 1, 12)


plotdata <- rbind(as.data.frame(t(rbind.data.frame(tMirExp,tGeneExp))) %>% mutate(Type = 'Tumor'),
as.data.frame(t(rbind.data.frame(nMirExp,nGeneExp))) %>% mutate(Type = 'Normal'))
colnames(plotdata) <- c('miRNA', 'Target', 'Type')



p0 <- ggplot(plotdata, aes(miRNA, Target, colour = Type)) +
  geom_point() + theme_bw() + 
  scale_color_manual(values = brewer.pal(3,'Set2')[1:2]) + 
  theme(legend.position = c(0.12,0.12), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.background = element_blank(), legend.title = element_blank()) + 
  labs(x = "Expression of hsa-miR-339-5p", y = "Expression of P2RY2")

miRNATargetPairedExpPlot <- ggMarginal(p0, type = "boxplot", groupColour = TRUE, groupFill = TRUE)


############
ggsave(plot = miRNATargetExpPlot + miRNATargetPairedExpPlot, width = 14, height = 7, units = 'cm',
       filename = 'D:/CancerNeuroscience/Github/result/section4/miR-339_P2RY2_PRAD_ExpPlot.pdf')

