

library('dplyr')
library('ggplot2')
library('ggExtra')
library('Hmisc')
library('RColorBrewer')
##########
tcgaPanCanSamples <- readRDS(file = 'D:/CancerNeuroscience/Github/data/tcgaPanCanSamples.rds')

load(file = 'D:/CancerNeuroscience/Github/data/panCanGeneExpData.RData')
# geneInfo, panCanTurGeneExp, panCanPairdTurGeneExp, panCanPairdNormGeneExp

load(file = 'D:/CancerNeuroscience/Github/data/panCanMethyData.RData')
# panCanTurMethy, panCanPairdTurMethy, panCanPairdNormMethy


##############
anno450k <- readRDS(file = 'D:/CancerNeuroscience/Github/data/methAnno450k.rds')
nrAnno450k <- subset(anno450k, RefGene_Name %in% 'GRIA1')

nrpromoterAnno <- subset(nrAnno450k, Name %in% c('cg17020834'))

############## Methy
panCanTurMethy <- rownames_to_column(panCanTurMethy, var = 'Name')
TMethy <- panCanTurMethy %>% inner_join(nrpromoterAnno[, c('Name', 'RefGene_Name')], by = 'Name')

TMethy <- TMethy[!duplicated(TMethy$Name), ]
TMethy <- TMethy %>%  select(-RefGene_Name) %>% remove_rownames() %>% column_to_rownames(var = 'Name')  
TMethy <- colMeans(TMethy) %>% as.data.frame() %>% t()%>% as.data.frame()

rownames(TMethy) <- 'Methy'


panCanPairdNormMethy <- rownames_to_column(panCanPairdNormMethy, var = 'Name')
NMethy <- panCanPairdNormMethy %>% inner_join(nrpromoterAnno[, c('Name', 'RefGene_Name')], by = 'Name')

NMethy <- NMethy[!duplicated(NMethy$Name), ]
NMethy <- NMethy %>%  select(-RefGene_Name) %>% remove_rownames() %>% column_to_rownames(var = 'Name')  
NMethy <- colMeans(NMethy) %>% as.data.frame() %>% t()%>% as.data.frame()

rownames(NMethy) <- 'Methy'

############## Exp
TExp <- panCanTurGeneExp['2890', ]
NExp <- panCanPairdNormGeneExp['2890', ]
rownames(TExp) <- 'Exp'
rownames(NExp) <- 'Exp'


##############
analySamples <- subset(tcgaPanCanSamples, DISEASE == 'LUAD')
comSmas <- intersect(intersect(colnames(TExp), colnames(TMethy)), analySamples$SAMPLE_BARCODE)
methyGeneExp <- as.data.frame(t(rbind.data.frame(TMethy[, comSmas], TExp[, comSmas])))


methyExpPlot <- ggplot(methyGeneExp, aes(x = Methy, y = Exp)) + geom_point() + geom_smooth(method = "glm", se = FALSE) + 
  labs(x = "Methylation of GRIA1 \n cg17020834", y = "Expression of GRIA1") + theme_bw() + 
  theme(legend.position = c(0.12,0.12), legend.background = element_blank(), legend.title = element_blank())


############
comSmas <- Reduce(intersect, list(substr(colnames(TMethy), 1, 12), substr(colnames(NMethy), 1, 12),
                                  substr(colnames(TExp), 1, 12), substr(colnames(NExp), 1, 12)))

sams <- subset(tcgaPanCanSamples, DISEASE == 'LUAD' & PATIENT_BARCODE %in% comSmas)



tGeneMethy <- TMethy[, paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tGeneMethy) <- substr(colnames(tGeneMethy), 1, 12)

tGeneExp <- TExp[, paste0(sams$PATIENT_BARCODE, '-01')]
colnames(tGeneExp) <- substr(colnames(tGeneExp), 1, 12)

nGeneMethy <- NMethy[, paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nGeneMethy) <- substr(colnames(nGeneMethy), 1, 12)

nGeneExp <- NExp[, paste0(sams$PATIENT_BARCODE, '-11')]
colnames(nGeneExp) <- substr(colnames(nGeneExp), 1, 12)


plotdata <- rbind(as.data.frame(t(rbind.data.frame(tGeneMethy,tGeneExp))) %>% mutate(Type = 'Tumor'),
                  as.data.frame(t(rbind.data.frame(nGeneMethy,nGeneExp))) %>% mutate(Type = 'Normal'))



p0 <- ggplot(plotdata, aes(Methy, Exp, colour = Type)) + geom_point() + theme_bw() + 
  scale_color_manual(values = brewer.pal(3,'Set2')[1:2]) + 
  labs(x = "Methylation of GRIA1 \n cg17020834", y = "Expression of GRIA1")+ 
  theme(legend.position = c(0.12,0.12), legend.background = element_blank(), legend.title = element_blank(),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7))

methyExpPairedPlot <- ggMarginal(p0, type = "boxplot", groupColour = TRUE, groupFill = TRUE)


############
ggsave(plot = methyExpPlot + methyExpPairedPlot, width = 14, height = 7, units = 'cm',
       filename = 'D:/CancerNeuroscience/Github/result/section4/cg17020834_GRIA1_LUAD_ExpPlot.pdf')

