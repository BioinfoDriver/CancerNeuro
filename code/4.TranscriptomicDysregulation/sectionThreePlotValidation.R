library(ggplot2)
library(patchwork)
library(pheatmap)
library(ggpubr)
library(dplyr)
library(scatterpie)
library(tibble)

#########################
nRDiffExp <- readRDS(file = '/data/curatedPanCanNrDiffExp.rds')


diseasedDysGStat <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% group_by(DISEASE) %>% 
  dplyr::count(name = 'numOfDysGenes')

diseaseUpGStat <- nRDiffExp %>% subset(logFC > 1 & adj.P.Val < 0.05) %>% group_by(DISEASE) %>% 
  dplyr::count() %>% mutate(DisType = 'Up')
diseaseDownGStat <- nRDiffExp %>% subset(logFC < -1 & adj.P.Val < 0.05) %>% group_by(DISEASE) %>% 
  dplyr::count() %>% mutate(DisType = 'Down') 


percOfDysGS <- rbind.data.frame(diseaseUpGStat, diseaseDownGStat) %>% 
  left_join(diseasedDysGStat, by = join_by(DISEASE)) %>% mutate(percentage = n/numOfDysGenes)



plot1 <- ggplot(percOfDysGS, aes(x = "", y = percentage, fill = DisType)) +
  geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) +
  theme_void() + facet_wrap(~ DISEASE, nrow = 2) + 
  scale_fill_manual(values = c('Down' = 'blue', 'Up' = "red"), name = '') + 
  geom_text(aes(label = paste0(round(percentage * 100, 0), "%")), 
            position = position_stack(vjust = 0.5), size = 3, colour = 'black')

#########################
# percOfDysGS$x_locus <- rep(c(seq(2, 22, 4), seq(2, 22, 4)), 2)
# percOfDysGS$y_locus <- rep(c(rep(2, 6), rep(6, 6)), 2)
# 
# percOfDysGS$value <- percOfDysGS$percentage
# percOfDysGS$DisType <- factor(percOfDysGS$DisType, levels = c('Down', 'Up'))
# 
# plot2 <- ggplot() + geom_scatterpie(aes(x=x_locus, y=y_locus, r=log10(numOfDysGenes* 0.5)), 
#                                     data=percOfDysGS, cols="DisType", long_format=TRUE) + coord_fixed() + 
#   scale_fill_manual(values = c('Down' = 'blue', 'Up' = "red"), name = '') + 
#   ggrepel::geom_text_repel(aes(x_locus, y_locus, label = n), data=percOfDysGS, size = 3)  + 
#   ggrepel::geom_text_repel(aes(x_locus, y_locus + 2, label = DISEASE), data=percOfDysGS, size = 3)


#########################

geneDysGStat <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% group_by(geneSymbol) %>% 
  dplyr::count(name = 'numOfDisDys')

geneUpGStat <- nRDiffExp %>% subset(logFC > 1 & adj.P.Val < 0.05) %>% group_by(geneSymbol) %>% 
  dplyr::count() %>% mutate(DisType = 'Up')
geneDownGStat <- nRDiffExp %>% subset(logFC < -1 & adj.P.Val < 0.05) %>% group_by(geneSymbol) %>% 
  dplyr::count() %>% mutate(DisType = 'Down') 



geneDysGS <- rbind.data.frame(geneUpGStat, geneDownGStat) %>% left_join(geneDysGStat, by = join_by(geneSymbol)) %>%
  arrange(desc(numOfDisDys), desc(n)) %>% mutate(DisType = factor(DisType, levels = c('Up', 'Down')))



plot3 <- ggplot(subset(geneDysGS, numOfDisDys > 3), aes(x = reorder(geneSymbol, numOfDisDys), y = n, fill = DisType)) +
  geom_bar(stat = "identity", position = "stack") + scale_y_continuous(breaks = seq(0, 12, 3), position = "right") +
  labs(y = 'Number of cancers') + 
  scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), name = '') +
  theme(axis.title.y = element_blank(), 
        legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black")) + coord_flip()


#########################

logFcHeatmap <- reshape2::dcast(nRDiffExp, geneSymbol ~ DISEASE, value.var = "logFC")
logFcHeatmap <- column_to_rownames(logFcHeatmap, var = "geneSymbol")
# logFcHeatmap[is.na(logFcHeatmap)] <- 0


nRDiffExp <- nRDiffExp %>% mutate(diffIndex = abs(logFC) > 1 & adj.P.Val < 0.05)
borderColors <- reshape2::dcast(nRDiffExp, geneSymbol ~ DISEASE, value.var = "diffIndex")

borderColors <- column_to_rownames(borderColors, var = "geneSymbol")
borderColors <- t(apply(borderColors, 1, function(x) ifelse(x == TRUE, 'black', NA)))


showGene <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% group_by(geneSymbol) %>% 
  dplyr::count(name = 'numOfDisDys') %>% arrange(desc(numOfDisDys)) %>% subset(numOfDisDys > 3, 'geneSymbol')
showGene <- showGene$geneSymbol


plot4 <- pheatmap(logFcHeatmap[showGene, ], cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
                  show_rownames = TRUE, number_color = "grey30", number_format = "%.2f", show_colnames = TRUE, 
                  angle_col = 45, color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_number = 6, fontsize = 7, 
                  cellheight = 12, cellwidth = 24, border = borderColors[showGene, ])


#########################


diseaseDysGStat <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% mutate(DisType = ifelse(logFC > 1, 'Up', 'Down')) %>% 
  group_by(DISEASE, DisType) %>% dplyr::count(name = 'numOfDisDys')


plot5 <- ggplot(diseaseDysGStat, aes(x = DISEASE, y = numOfDisDys, fill = DisType)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(y = "Number of dysregulated gene") +
  scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), name = '') +
  theme(legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.title.y = element_blank(), axis.text.x = element_text(angle = -270, hjust = 0.5, vjust = 0.5)) + coord_flip()


#########################################
pdf(file = '/result/section3/validation/Figure3.1.pdf')

plot1
plot3 + plot5

dev.off()

pdf(file = '/result/section3/validation/Figure3.2.pdf')

plot4

dev.off()


