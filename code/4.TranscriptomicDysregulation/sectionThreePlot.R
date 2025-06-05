
library(ggplot2)
library(patchwork)
library(pheatmap)
library(ggpubr)
library(dplyr)
library(scatterpie)
library(tibble)

#########################
nRDiffExp <- readRDS(file = '/data/panCanNrDiffExp.rds')


diseasedDysGStat <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% group_by(DISEASE) %>% 
  count(name = 'numOfDysGenes')

diseaseUpGStat <- nRDiffExp %>% subset(logFC > 1 & adj.P.Val < 0.05) %>% group_by(DISEASE) %>% 
  count() %>% mutate(DisType = 'Up')
diseaseDownGStat <- nRDiffExp %>% subset(logFC < -1 & adj.P.Val < 0.05) %>% group_by(DISEASE) %>% 
  count() %>% mutate(DisType = 'Down') 


percOfDysGS <- rbind.data.frame(diseaseUpGStat, diseaseDownGStat) %>% 
  left_join(diseasedDysGStat, by = join_by(DISEASE)) %>% mutate(percentage = n/numOfDysGenes)



plot1 <- ggplot(percOfDysGS, aes(x = "", y = percentage, fill = DisType)) +
  geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) +
  theme_void() + facet_wrap(~ DISEASE, nrow = 2) + 
  scale_fill_manual(values = c('Down' = 'blue', 'Up' = "red"), name = '') + 
  geom_text(aes(label = paste0(round(percentage * 100, 0), "%")), 
            position = position_stack(vjust = 0.5), size = 3, colour = 'black')

#########################
percOfDysGS$x_locus <- rep(c(seq(2, 24, 3), seq(2, 21, 3)), 2)
percOfDysGS$y_locus <- rep(c(rep(2, 8), rep(6, 7)), 2)

percOfDysGS$value <- percOfDysGS$percentage
percOfDysGS$DisType <- factor(percOfDysGS$DisType, levels = c('Down', 'Up'))

plot2 <- ggplot() + geom_scatterpie(aes(x=x_locus, y=y_locus, r=log10(numOfDysGenes* 0.5)), 
                                    data=percOfDysGS, cols="DisType", long_format=TRUE) + coord_fixed() + 
  scale_fill_manual(values = c('Down' = 'blue', 'Up' = "red"), name = '') + 
  ggrepel::geom_text_repel(aes(x_locus, y_locus, label = n), data=percOfDysGS, size = 3)  + 
  ggrepel::geom_text_repel(aes(x_locus, y_locus + 2, label = DISEASE), data=percOfDysGS, size = 3)


#########################

geneDysGStat <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% group_by(Approved.symbol) %>% 
  count(name = 'numOfDisDys')

geneUpGStat <- nRDiffExp %>% subset(logFC > 1 & adj.P.Val < 0.05) %>% group_by(Approved.symbol) %>% 
  count() %>% mutate(DisType = 'Up')
geneDownGStat <- nRDiffExp %>% subset(logFC < -1 & adj.P.Val < 0.05) %>% group_by(Approved.symbol) %>% 
  count() %>% mutate(DisType = 'Down') 



geneDysGS <- rbind.data.frame(geneUpGStat, geneDownGStat) %>% left_join(geneDysGStat, by = join_by(Approved.symbol)) %>%
  arrange(desc(numOfDisDys), desc(n)) %>% mutate(DisType = factor(DisType, levels = c('Up', 'Down')))



plot3 <- ggplot(subset(geneDysGS, numOfDisDys > 7), aes(x = reorder(Approved.symbol, numOfDisDys), y = n, fill = DisType)) +
  geom_bar(stat = "identity", position = "stack") + scale_y_continuous(breaks = seq(0, 12, 3), position = "right") +
  labs(y = 'Number of cancers') + 
  scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), name = '') +
  theme(axis.title.y = element_blank(), 
        legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black")) + coord_flip()


#########################

logFcHeatmap <- reshape2::dcast(nRDiffExp, Approved.symbol ~ DISEASE, value.var = "logFC")
logFcHeatmap <- column_to_rownames(logFcHeatmap, var = "Approved.symbol")
# logFcHeatmap[is.na(logFcHeatmap)] <- 0


nRDiffExp <- nRDiffExp %>% mutate(diffIndex = abs(logFC) > 1 & adj.P.Val < 0.05)
borderColors <- reshape2::dcast(nRDiffExp, Approved.symbol ~ DISEASE, value.var = "diffIndex")

borderColors <- column_to_rownames(borderColors, var = "Approved.symbol")
borderColors <- t(apply(borderColors, 1, function(x) ifelse(x == TRUE, 'black', NA)))


showGene <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% group_by(Approved.symbol) %>% 
  count(name = 'numOfDisDys') %>% arrange(desc(numOfDisDys)) %>% subset(numOfDisDys > 7, 'Approved.symbol')
showGene <- showGene$Approved.symbol


plot4 <- pheatmap(logFcHeatmap[showGene, ], cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
                  show_rownames = TRUE, number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, 
                  angle_col = 45, color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_number = 6, fontsize = 7, 
                  cellheight = 12, cellwidth = 24, border = borderColors[showGene, ])


#########################


diseaseDysGStat <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% mutate(DisType = ifelse(logFC > 1, 'Up', 'Down')) %>% 
  group_by(DISEASE, DisType) %>% count(name = 'numOfDisDys')


plot5 <- ggplot(diseaseDysGStat, aes(x = DISEASE, y = numOfDisDys, fill = DisType)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(y = "Number of dysregulated gene") +
  scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), name = '') +
  theme(legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.title.y = element_blank(), axis.text.x = element_text(angle = -270, hjust = 0.5, vjust = 0.5)) + coord_flip()
# 
# 
# pdf(file = '/result/section3/new/Figure3.1.pdf')
# 
# plot1 / plot2
# plot3 + plot5
# 
# dev.off()
# 
# pdf(file = '/result/section3/new/Figure3.2.pdf')
# 
# plot4
# 
# dev.off()
# 


#########################
load(file = '/data/panCanGeneExpData.RData')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

panCanPairdNormGeneExp <- panCanPairdNormGeneExp[as.character(neurotransmitterReceptors$NCBI.Gene.ID), ]
panCanPairdTurGeneExp <- panCanPairdTurGeneExp[as.character(neurotransmitterReceptors$NCBI.Gene.ID), ]

rownames(panCanPairdNormGeneExp) <- neurotransmitterReceptors$Approved.symbol
rownames(panCanPairdTurGeneExp) <- neurotransmitterReceptors$Approved.symbol

panCanPairdNormGeneExp <- as.data.frame(t(panCanPairdNormGeneExp))
panCanPairdTurGeneExp <- as.data.frame(t(panCanPairdTurGeneExp))


rownames(panCanPairdNormGeneExp) <- substr(rownames(panCanPairdNormGeneExp), 1, 12)
rownames(panCanPairdTurGeneExp) <- substr(rownames(panCanPairdTurGeneExp), 1, 12)
panCanPairdTurGeneExp <- panCanPairdTurGeneExp[rownames(panCanPairdNormGeneExp), ]


########
showData <- data.frame(normalExp = panCanPairdNormGeneExp[, 'ADRA1A'], tumorExp = panCanPairdTurGeneExp[, 'ADRA1A'],
                       PATIENT_BARCODE = rownames(panCanPairdNormGeneExp)) %>% 
  left_join(tcgaPanCanSamples, by = join_by(PATIENT_BARCODE))


showData <- showData %>% group_by(DISEASE) %>% count() %>% subset(n >=10) %>% left_join(showData, by = join_by(DISEASE))

colPalette <- c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
             '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3', '#D92E88', '#4D2A80', '#B2202B', '#D07927', '#CFBFDC', '#0B477A', 
             '#791A1C', '#A64C93', '#967FB4', '#0FA095', '#1EA2DC', '#DEBC27', '#EEE33E', '#1177A9', '#69779A', '#EEAAAE', '#CB96BF', '#ED8F27')

plot6.1 <- ggpaired(showData, cond1 = "normalExp", cond2 = "tumorExp", facet.by = 'DISEASE', fill = "DISEASE", color = "black",
         line.color = "gray90", line.size = 0.4, palette = colPalette, xlab = FALSE, 
         ylab = 'Normalized Gene Expression')



showData <- data.frame(normalExp = panCanPairdNormGeneExp[, 'GABRD'], tumorExp = panCanPairdTurGeneExp[, 'GABRD'],
                       PATIENT_BARCODE = rownames(panCanPairdNormGeneExp)) %>% 
  left_join(tcgaPanCanSamples, by = join_by(PATIENT_BARCODE))

showData <- showData %>% group_by(DISEASE) %>% count() %>% subset(n >=10) %>% left_join(showData, by = join_by(DISEASE))

plot6.2 <- ggpaired(showData, cond1 = "normalExp", cond2 = "tumorExp", facet.by = 'DISEASE', fill = "DISEASE", color = "black",
                    line.color = "gray90", line.size = 0.4, palette = colPalette, xlab = FALSE, 
                    ylab = 'Normalized Gene Expression')


showData <- data.frame(normalExp = panCanPairdNormGeneExp[, 'CHRNA5'], tumorExp = panCanPairdTurGeneExp[, 'CHRNA5'],
                       PATIENT_BARCODE = rownames(panCanPairdNormGeneExp)) %>% 
  left_join(tcgaPanCanSamples, by = join_by(PATIENT_BARCODE))

showData <- showData %>% group_by(DISEASE) %>% count() %>% subset(n >=10) %>% left_join(showData, by = join_by(DISEASE))

plot6.3 <- ggpaired(showData, cond1 = "normalExp", cond2 = "tumorExp", facet.by = 'DISEASE', fill = "DISEASE", color = "black",
                    line.color = "gray90", line.size = 0.4, palette = colPalette, xlab = FALSE, 
                    ylab = 'Normalized Gene Expression')


showData <- data.frame(normalExp = panCanPairdNormGeneExp[, 'HTR1D'], tumorExp = panCanPairdTurGeneExp[, 'HTR1D'],
                       PATIENT_BARCODE = rownames(panCanPairdNormGeneExp)) %>% 
  left_join(tcgaPanCanSamples, by = join_by(PATIENT_BARCODE))

showData <- showData %>% group_by(DISEASE) %>% count() %>% subset(n >=10) %>% left_join(showData, by = join_by(DISEASE))

plot6.4 <- ggpaired(showData, cond1 = "normalExp", cond2 = "tumorExp", facet.by = 'DISEASE', fill = "DISEASE", color = "black",
                    line.color = "gray90", line.size = 0.4, palette = colPalette, xlab = FALSE, 
                    ylab = 'Normalized Gene Expression')


# 
# pdf(file = '/result/section3/Figure3.3.pdf')
# 
# plot6.1
# plot6.2
# plot6.3
# plot6.4
# 
# dev.off()
# 

#########################

pairdSams <- rownames(panCanPairdTurGeneExp)
pairdSams <- tcgaPanCanSamples %>% subset(PATIENT_BARCODE %in% pairdSams)


samSta <- pairdSams %>% group_by(DISEASE) %>% count() %>% subset(n >=10)
pairdSams <- subset(pairdSams, DISEASE %in% samSta$DISEASE) %>% group_by(DISEASE) %>% count(name = 'numOfPiredSams') 


fillCol <- c('#734B27', '#364D99', '#168C42', '#B9A131', '#C7CAD9', '#C9A88D', '#EDB067', '#B0CC46', '#F5CED6', 
             '#DE6B71', '#98D6F0', '#92C9A4', '#D7EDF7', '#DC2425', '#F6DFC3')


plot7 <- ggplot(pairdSams, aes(x = reorder(DISEASE, -numOfPiredSams), y = numOfPiredSams)) +
  geom_bar(stat = "identity", fill = fillCol)  + labs(y = 'Number of Patients') + geom_text(aes(label = numOfPiredSams), vjust = -0.5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 


# 
# pdf(file = '/result/section3/SFigure3.1.pdf')
# 
# plot7
# 
# dev.off()
# 


