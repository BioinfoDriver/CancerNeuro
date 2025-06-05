
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)

load(file = '/data/neuroRecepCnStat.RData')
load(file = '/data/neuroRecepCnComToRandom.RData')

# nRCnStatByDis, nRCnStatByGroupDis, nRCnStatByGeneDis, nRCnStatByGene,

##############

nrCnaGain <- subset(nRCnStatByGroupDis, cnStaus == 2) %>% 
  reshape2::dcast(classesOfNeuroreceptors ~ DISEASE, value.var = 'altFre', fill = 0) %>% 
  column_to_rownames(var = 'classesOfNeuroreceptors')

nrCnaLoss <- subset(nRCnStatByGroupDis, cnStaus == -2) %>% 
  reshape2::dcast(classesOfNeuroreceptors ~ DISEASE, value.var = 'altFre', fill = 0) %>% 
  column_to_rownames(var = 'classesOfNeuroreceptors')

nrCnaLoss$KICH <- 0
nrCnaLoss <- nrCnaLoss[, colnames(nrCnaGain)]

###
nrCnaGainPvalue <- comResByDisGroupPvalue[grep('_2', rownames(comResByDisGroupPvalue)), ]
rownames(nrCnaGainPvalue) <- gsub('_2', '', rownames(nrCnaGainPvalue))
nrCnaGainPvalue <- as.data.frame(t(nrCnaGainPvalue))

nrCnaLossPvalue <- comResByDisGroupPvalue[grep('_-2', rownames(comResByDisGroupPvalue)), ]
rownames(nrCnaLossPvalue) <- gsub('_-2', '', rownames(nrCnaLossPvalue))
nrCnaLossPvalue <- as.data.frame(t(nrCnaLossPvalue))

###

UpColor <- colorRamp2(breaks = c(0, 1), colors = c("white","#AB221F"))
DnColor <- colorRamp2(breaks = c(0, 1), colors = c("white","#3878C1"))


DiagFunc <- function(up, down, lossp, gainp){
  function(j, i, x, y, width, height, fill){

    
    if(lossp[i, j] < 0.05){
      
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "black"))
      
      grid.text('*', x - 0.25*width, y + 0.25*height, gp = gpar(fontsize = 15))
      
    }else{
      
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey"))
      
    }
    
    
    if(gainp[i, j] < 0.05){
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "black"))
      
      grid.text('*', x + 0.25*width, y - 0.25*height, gp = gpar(fontsize = 15))
      
    }else{
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    }

  }
}

p1 <- Heatmap(nrCnaGain, column_title = "Copy number variation across cancer types",
              rect_gp = gpar(type = "none"),
              show_heatmap_legend = F,
              cluster_rows = F,
              cluster_columns = F, 
              cell_fun = DiagFunc(up = nrCnaGain, down = nrCnaLoss, lossp = nrCnaLossPvalue, gainp = nrCnaGainPvalue))



lgd <- list(Legend(title = "CNV Gain", 
                   col_fun = UpColor, 
                   at = c(0,0.25, 0.5, 0.75, 1), 
                   direction = "horizontal" 
),
Legend(title = "CNV Loss", 
       col_fun = DnColor, 
       at = c(0,0.25, 0.5, 0.75, 1), 
       direction = "horizontal" 
))


pdf(file = '/result/section2/nrGroupCnaFre.pdf')
draw(p1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)

dev.off()


##############

comResByDisPvalue <- comResByDisPvalue %>% mutate(cnStaus = ifelse(cnStaus == 2, 'Amplification', 'Deletion')) %>% 
  mutate(altFre = ifelse(cnStaus == 'Amplification', altFre, -altFre))

sortIndex <- comResByDisPvalue %>% reshape2::dcast(DISEASE ~ cnStaus, value.var = "altFre", fill = 0) %>% 
  arrange(desc(Amplification)) %>% mutate(index = seq(length(Amplification)))


comResByDisPvalue <- comResByDisPvalue %>% left_join(sortIndex, by = 'DISEASE')



pdf(file = '/result/section2/nrTotalCnaFre.pdf')
ggplot(comResByDisPvalue, aes(x = DISEASE, y = altFre, fill = cnStaus)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") + theme_minimal() + 
  scale_fill_manual(values = c('Deletion' = "#3878C1", 'Amplification' = "#AB221F"), name = '') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "top", legend.justification = "left") + 
  geom_text(data = comResByDisPvalue,
            size = 3, check_overlap = TRUE, aes(label = pValue), nudge_y = 0.02)

dev.off()


##############
nRCnStatByGeneDis <- nRCnStatByGeneDis %>% mutate(cnStaus = ifelse(cnStaus == 2, 'Amplification', 'Deletion')) 

showGene <- unique(subset(nRCnStatByGeneDis, altFre>0.07)$Approved.symbol)


nrCnaGain <- subset(nRCnStatByGeneDis, cnStaus == 'Amplification') %>% 
  reshape2::dcast(Approved.symbol ~ DISEASE, value.var = 'altFre', fill = 0) %>% 
  column_to_rownames(var = 'Approved.symbol')

nrCnaLoss <- subset(nRCnStatByGeneDis, cnStaus == 'Deletion') %>% 
  reshape2::dcast(Approved.symbol ~ DISEASE, value.var = 'altFre', fill = 0) %>% 
  column_to_rownames(var = 'Approved.symbol')


nrCnaLoss$KICH <- 0
nrCnaLoss <- nrCnaLoss[, colnames(nrCnaGain)]


nrCnaGain <- nrCnaGain[showGene, ]
nrCnaLoss <- nrCnaLoss[showGene, ]



UpColor <- colorRamp2(breaks = c(0, 1), colors = c("white","#AB221F"))
DnColor <- colorRamp2(breaks = c(0, 1), colors = c("white","#3878C1"))


DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    
   if(down[i, j] >= 0.08){
     
     grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                  unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                  gp = gpar(fill = DnColor(down[i, j]), col = "grey"))
     
     grid.text(sprintf("%.2f", down[i, j]), x - 0.25*width, y + 0.25*height, gp = gpar(fontsize = 10))
     
   }else{
     
     grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                  unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                  gp = gpar(fill = DnColor(down[i, j]), col = "grey"))    
     
   } 


    if(up[i, j] >= 0.08){
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))      
      
      grid.text(sprintf("%.2f", up[i, j]), x + 0.25*width, y - 0.25*height, gp = gpar(fontsize = 10))
      
    }else{
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))         
      
    }  
      
      
  }
}

p1 <- Heatmap(nrCnaGain, column_title = "Copy number variation across cancer types",
              rect_gp = gpar(type = "none"),
              show_heatmap_legend = F,
              cluster_rows = F,
              cluster_columns = F, 
              cell_fun = DiagFunc(up = nrCnaGain, down = nrCnaLoss))



lgd <- list(Legend(title = "CNV Gain", 
                   col_fun = UpColor, 
                   at = c(0,0.25, 0.5, 0.75, 1), 
                   direction = "horizontal" 
),
Legend(title = "CNV Loss", 
       col_fun = DnColor, 
       at = c(0,0.25, 0.5, 0.75, 1), 
       direction = "horizontal" 
))

pdf(file = '/result/section2/nrGeneCnaFre.pdf')
draw(p1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)

dev.off()

##############
nRCnStatByGene <- nRCnStatByGene %>% mutate(cnStaus = ifelse(cnStaus == 2, 'Amplification', 'Deletion')) %>% 
  mutate(altFre = ifelse(cnStaus == 'Amplification', altFre, -altFre))

# showGene <- subset(nRCnStatByGene, abs(altFre) > 0.009)$Approved.symbol


tmp <- subset(nRCnStatByGene, Approved.symbol %in% showGene)
tmp$Approved.symbol <- factor(tmp$Approved.symbol, levels = rownames(nrCnaGain))


pdf(file = '/result/section2/nrGeneTotalCnaFre.pdf')

ggplot(tmp, 
       aes(x = Approved.symbol, y = altFre, fill = cnStaus)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") + theme_minimal() + 
  scale_fill_manual(values = c('Deletion' = "#3878C1", 'Amplification' = "#AB221F"), name = '') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "top", legend.justification = "left")

dev.off()


