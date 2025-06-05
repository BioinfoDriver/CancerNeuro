
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggbreak)

#########################

geneInfo <- readRDS(file = '/data/geneInfo.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')

neurotransmitterReceptors <- merge(neurotransmitterReceptors, 
                                   select(geneInfo, gene_name, exon_breaks), by.x = 'Approved.symbol', by.y = 'gene_name')

groupGeneNum <- neurotransmitterReceptors %>% group_by(classesOfNeuroreceptors) %>% count(name = 'NumOfGenes') %>%
  arrange(desc(NumOfGenes))
geneLengthstat <- neurotransmitterReceptors %>% group_by(classesOfNeuroreceptors, exon_breaks) %>% count(name = 'NumOfGenes')


geneLengthstat <- merge(geneLengthstat, groupGeneNum, by = 'classesOfNeuroreceptors')
geneLengthstat <- geneLengthstat %>% mutate(Percentage = NumOfGenes.x/NumOfGenes.y)
geneLengthstat$exon_breaks <- droplevels(geneLengthstat$exon_breaks)


##########
fillCol <- c('#F5E09B', '#DFE1E2', '#B7DBE3', '#9EE092', '#F5D9E6', '#D2D2D2', '#F4DEBB', '#E1B6B5', '#F9F2C1')
names(fillCol) <- groupGeneNum$classesOfNeuroreceptors
plot1 <- ggplot(groupGeneNum, aes(x = reorder(classesOfNeuroreceptors, -NumOfGenes), y = NumOfGenes)) +
  geom_bar(stat = "identity", fill = fillCol)  + labs(y = 'Number of genes') + 
  scale_fill_manual(values = fillCol) + geom_text(aes(label = NumOfGenes), vjust = -0.5) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  labs(tag ="(A)", color = "Species")

##########
plot2 <- ggplot(geneLengthstat, aes(x = reorder(classesOfNeuroreceptors, -NumOfGenes.y), y = Percentage, fill = exon_breaks)) +
  geom_bar(stat = "identity", position = "stack", color = "gray60")  + labs(y = 'Percentage of genes') + 
  scale_fill_manual(values = colorRampPalette(c("navy", "white", "firebrick3"))(36), name = 'Length') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", legend.key.size = unit(0.3, "cm")) + 
  labs(tag ="(B)", color = "Species")

##########
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

disCount <- mc3MutData %>% group_by(DISEASE, SAMPLE_BARCODE) %>% count(name = 'samMutNum') %>% 
  group_by(DISEASE) %>% count(name = 'numOfPat')


fillCol <- c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
             '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', '#D07927', '#CFBFDC', '#0B477A', 
             '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', '#EEE33E', '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27')

plot3 <- ggplot(disCount, aes(x = reorder(DISEASE, -numOfPat), y = numOfPat)) +
  geom_bar(stat = "identity", fill = fillCol)  + labs(y = 'Number of patients') + 
  geom_text(aes(label = numOfPat), vjust = -0.5, size = 4) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  labs(tag ="(C)", color = "Species")


##########
# plot4 <- plot1+plot2
# plot5 <- plot4 / plot3 + 
#   patchwork::plot_layout(design = "
#                          AAAA
#                          cccc
#                          ", widths = c(1, 1))
# ggsave(filename = '/result/section1/SFigure1.pdf', plot = plot5)
# 

# pdf('/result/section1/SFigure1.pdf')
# 
# plot1
# plot2
# plot3
# 
# dev.off()


#########################
gainLoFMutStat <- readRDS(file = '/data/gainLoFMutStat.rds')

# size
gainLoFMutStat$sizeLable <- cut(gainLoFMutStat$numOfNonsMut, breaks = c(0, 32, 50, 86, 300), 
                                labels = c('n ≤ 32', '32 < n ≤ 50', '50 < n ≤ 86', 'n > 86'))

# LoF OR Hotspot
gainLoFMutStat$colLabel <- ifelse(gainLoFMutStat$gainMf > 0.3 & gainLoFMutStat$lofMf < 0.2 & gainLoFMutStat$numOfHotMutPosi >5, 'HotspotMutG', 
       ifelse(gainLoFMutStat$gainMf < 0.2 & gainLoFMutStat$lofMf > 0.3 & gainLoFMutStat$numOfLofMut > 5, 'LofMutG', 'OtherG'))


gainLoFMutStat <- tibble::rownames_to_column(gainLoFMutStat, var = "GeneName")

colLabel <- c("#B7DBE3", "#DFE1E2", "#F5D9E6")
names(colLabel) <- c('LofMutG', 'OtherG', 'HotspotMutG')

sizeLable <- c(1.5, 3, 4.5, 6)
names(sizeLable) <- c('n ≤ 32', '32 < n ≤ 50', '50 < n ≤ 86', 'n > 86')


plot6 <- ggplot(gainLoFMutStat, aes(x = lofMf, y = gainMf, color = colLabel, size = sizeLable)) +
  geom_point() + labs(x = "Fraction of LoF mutations", y = "Fracton of Hotspot mutations") + 
  scale_colour_manual(values = colLabel, guide = "legend", name = "Color") + 
  scale_size_manual(values = sizeLable, guide = "legend", name = "Size") + 
  theme(legend.position = "top", legend.justification = "left", panel.background = element_blank()) + 
  geom_text(data = subset(gainLoFMutStat, colLabel %in% c('LofMutG', 'HotspotMutG')),
            size = 3, check_overlap = TRUE, aes(label = GeneName), nudge_y = 0.005) + 
  geom_hline(yintercept = -0.01, color = "black", linetype = "solid") + 
  geom_vline(xintercept = -0.01, color = "black", linetype = "solid") + 
  scale_x_break(breaks=c(0.3, 0.5), ticklabels = c(0.6))


##########################

load(file = '/data/neuroRecepMutStat.RData')
load(file = '/data/neuroRecepMutSig2CV.RData')


# install.packages("reshape2")
library(reshape2)

mutSigGenePvalue <- melt(mutSigGenePvalue, id.vars = 'gene', value.name = 'pValue', variable.name = 'cancerType')
mutSigGeneQvalue <- melt(mutSigGeneQvalue, id.vars = 'gene', value.name = 'qValue', variable.name = 'cancerType')

mutSigGenePvalue <- mutSigGenePvalue %>% mutate(cancerType = gsub('pvalue', '', cancerType))
mutSigGeneQvalue <- mutSigGeneQvalue %>% mutate(cancerType = gsub('qvalue', '', cancerType))


mutSigGenePvalue <- subset(mutSigGenePvalue, pValue < 0.01)

sigMut <- merge(mutSigGenePvalue, mutSigGeneQvalue, by = c('gene', 'cancerType'))

neuroRecepMutStatByGene$DISEASE = 'PANCAN'
neuroRecepMutStatByGene$numOfPat = 8217
neuroRecepMutStatByGene <- rename(neuroRecepMutStatByGene, numOfPatWithMut = numOfMut)
neuroRecepMutStatByGene <- neuroRecepMutStatByGene[, colnames(neuroRecepMutStatByDisGene)]

nrMutSat <- rbind.data.frame(neuroRecepMutStatByGene, neuroRecepMutStatByDisGene)

sigMut <- rename(sigMut, Hugo_Symbol = gene, DISEASE = cancerType)
sigMut <- merge(sigMut, nrMutSat, by = c('Hugo_Symbol', 'DISEASE'))



# size
sigMut$sizeLable <- cut(sigMut$pValue, breaks = c(1.0e-05, 1.0e-04, 1.0e-03, 1.0e-02), 
 labels = c('p < 1.0e-04','p < 1.0e-03', 'p < 1.0e-02'))


# shape
# sigMut$shapeLable <- cut(sigMut$qValue, breaks = c(0, 0.05, 0.25, 1.0), 
#                         labels = c('q < 0.05','q < 0.25', 'p < 1.0'))

sizeLable <- c(4, 5, 7)
names(sizeLable) <- c('p < 1.0e-02', 'p < 1.0e-03', 'p < 1.0e-04')

# shapeLable <- c( 17, 15, 16)
# names(shapeLable) <- c('q < 0.05','q < 0.25', 'p < 1.0')


# plot7 <- ggplot(sigMut, aes(x = Hugo_Symbol, y = DISEASE, color = altFre, size = sizeLable)) +
#   geom_point(aes(shape = shapeLable)) +
#   labs(x = "Significantly mutated genes", y = "Cancer Types") + 
#   scale_size_manual(values = sizeLable, guide = "legend", name = "MutSigCV p value") + 
#   scale_color_gradient(low = '#ECC38C', high = '#BD4146', name = 'Frequency') + 
#   scale_shape_manual(values = shapeLable, guide = "legend", name = "MutSigCV q value") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 


sigMut$colorLable <- cut(sigMut$qValue, breaks = c(0, 0.25, 1.0), 
                         labels = c('q < 0.25', 'q ≥ 0.25'))

colLabel <- c("black", "white")
names(colLabel) <- c('q < 0.25', 'q ≥ 0.25')


plot7 <- ggplot(sigMut, aes(x = Hugo_Symbol, y = DISEASE, fill = altFre, size = sizeLable, color = colorLable)) +
  geom_point(shape = 21, stroke = 1.5) +
  labs(x = "Significantly mutated genes", y = "Cancer Types") + 
  scale_size_manual(values = sizeLable, guide = "legend", name = "MutSigCV p value") + 
  scale_fill_gradient(low = '#ECC38C', high = '#BD4146', name = 'Frequency') + 
  scale_colour_manual(values = colLabel, guide = "legend", name = "MutSigCV q value") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(plot7, filename = '/result/section1/Fig1B_new.pdf')


##########################

load(file = '/data/neuroRecepMutStat.RData')
load(file = '/data/neuroRecepMutSig2CV.RData')
load('/data/neuroRecepMutComToRandomExonLength_Threshold_LargeThan.RData')
NonsilentTBM <- readRDS(file = '/data/meanOfNonsilentMutBurden.rds')


library(tibble)
library(pheatmap)

heatmapData <- reshape2::dcast(neuroRecepMutStatByDisGroup, DISEASE ~ classesOfNeuroreceptors, value.var = "altFre")
heatmapData <- column_to_rownames(heatmapData, var = "DISEASE")

heatmapData[is.na(heatmapData)] <- 0
heatmapData <- as.data.frame(t(heatmapData))


borderColors <- apply(comResByDisGroupPvalue, 1, function(x) ifelse(x < 0.05, 'black', NA))

NonsilentTBM <- NonsilentTBM %>% column_to_rownames(var = "DISEASE") %>% mutate(TBM = round(n, 1), .keep = 'unused')


plot8 <- pheatmap(heatmapData, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, show_rownames = TRUE, 
         number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, angle_col = 45, 
         fontsize_number = 5, fontsize = 7, 
         color = colorRampPalette(c("white", "firebrick3"))(50), cellheight = 12, cellwidth = 14, border = borderColors, 
         annotation_col = NonsilentTBM, annotation_colors = list(TBM = colorRampPalette(c("white", "#9EE092"))(100)))


##########################

nRMutStat <- neuroRecepMutStatByGene %>% mutate(DISEASE = 'ZPanCan', numOfPatWithMut = numOfMut, numOfPat = 8217)
nRMutStat <- rbind.data.frame(nRMutStat[, colnames(neuroRecepMutStatByDisGene)], neuroRecepMutStatByDisGene)


nRMutheatmap <- reshape2::dcast(nRMutStat, Hugo_Symbol~DISEASE, value.var = "altFre")
nRMutheatmap <- column_to_rownames(nRMutheatmap, var = "Hugo_Symbol")


nRMutheatmap[is.na(nRMutheatmap)] <- 0
nRMutheatmap <- nRMutheatmap[apply(nRMutheatmap, 1, function(x) sum(x > 0.01)>10), ]

borderColors <- t(apply(nRMutheatmap, 1, function(x) ifelse(x > 0.03, 'black', NA)))



plot9 <- pheatmap(nRMutheatmap, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, show_rownames = TRUE, 
  number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, angle_col = 45, 
  fontsize_number = 5, fontsize = 7, 
  color = colorRampPalette(c("white", "firebrick3"))(50), cellheight = 12, cellwidth = 15, border = borderColors)


# pdf(file = '/result/section1/SFigure2.pdf')
# plot8
# dev.off()
# 
# pdf(file = '/result/section1/ExonLengthSFigure2.pdf')
# plot8
# dev.off()

# pdf(file = '/result/section1/SFigure3.pdf')
# plot9
# dev.off()
# 

##########################

nRMutPvalue <- as.data.frame(comResByDisPvalue) %>% rownames_to_column(var = 'DISEASE') %>% 
  mutate(freq = 1 - comResByDisPvalue, .keep = 'unused') %>% arrange(desc(freq))

fillCol <- c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
  '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', '#D07927', '#CFBFDC', '#0B477A', 
  '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', '#EEE33E', '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27')

plot10 <- ggplot(nRMutPvalue, aes(x = reorder(DISEASE, -freq), y = freq)) +
  geom_bar(stat = "identity", fill = fillCol)  + labs(y = 'Normalized mutation load') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  geom_hline(yintercept = c(0.05, 0.95), color = c('black', "red"), linetype = "dashed", linewidth = 1) 


##########################

# plot11 <- plot6+plot7+plot10+
#   patchwork::plot_layout(design = "
#                          AABB
#                          AABB
#                          ####
#                          cccc")
# 
# 
# ggsave(filename = '/result/section1/Figure1.pdf', plot = plot11)
# 

# 
# pdf('/result/section1/Figure1.pdf')
# plot6
# plot7
# plot10
# dev.off()
# 


# 
# pdf('/result/section1/nRMutPvalueByDisease.pdf')
# plot10
# dev.off()


# pdf('/section1/ExonLength_Threshold_LargeThan.pdf')
# print(plot9)
# print(plot10)
# print(plot8)
# dev.off()


