
library('dplyr')
library('tibble')
library('ggplot2')
########################
miRnaDiffExp <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanMiRnaDiffExp.rds')
nRDiffExp <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanNrDiffExp.rds')
miRnaTarByDisease <- readRDS(file = 'D:/CancerNeuroscience/Github/data/nROfmiRNARegulatorFLM.rds')
panCanSurPvalue <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanRhoSurUniMulLogPvalue_PFI.rds')


#############
colnames(nRDiffExp) <- c('NR_geneID', 'NR_logFC', 'NR_AveExpr', 'NR_t', 'NR_P.Value', 'NR_adj.P.Val', 'NR_B', 'DISEASE', 'Symbol')
colnames(miRnaDiffExp) <- c('miRNA', 'miR_logFC', 'miR_AveExpr', 'miR_t', 'miR_P.Value', 'miR_adj.P.Val', 'miR_B', 'DISEASE')
colnames(miRnaTarByDisease) <- c("miRNA", "Estimate", "StdError", "tValue", "pValue", "DISEASE", "Symbol")
colnames(panCanSurPvalue)[1:2] <- c("DISEASE", "Symbol")


miRnaDiffExp <- miRnaDiffExp %>% mutate(MirStatus = ifelse(miR_logFC > 0.585 & miR_adj.P.Val < 0.05, 'Up', ifelse(miR_logFC < -0.585 & miR_adj.P.Val < 0.05, 'Down', 'Neutral')))
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(NR_logFC > 1 & NR_adj.P.Val < 0.05, 'Up', ifelse(NR_logFC < -1 & NR_adj.P.Val < 0.05, 'Down', 'Neutral')))
miRnaTarByDisease <- subset(miRnaTarByDisease, pValue < 0.05 & Estimate < -0.15)
panCanSurPvalue <- panCanSurPvalue %>% subset(mulcoxQvalue <= 0.1)


miRnaDiffExp <- miRnaDiffExp %>% right_join(miRnaTarByDisease, by = c('DISEASE', 'miRNA'))
nRDiffExp <- nRDiffExp %>% left_join(miRnaDiffExp,  by = join_by(DISEASE, Symbol)) 
nRDiffExp <- subset(nRDiffExp, !is.na(MirStatus)) %>% 
  subset((Estatus == 'Up' & MirStatus == 'Down') | (Estatus == 'Down' & MirStatus == 'Up'))

nRDiffExp <- nRDiffExp %>% left_join(panCanSurPvalue,  by = join_by(DISEASE, Symbol)) 

# na.omit(nRDiffExp)
#    NR_geneID  NR_logFC NR_AveExpr      NR_t   NR_P.Value NR_adj.P.Val      NR_B DISEASE Symbol Estatus           miRNA miR_logFC miR_AveExpr    miR_t
# 64      2903 -1.030606   4.012019 -4.826753 1.330815e-05 5.185887e-05  2.683808    PRAD GRIN2A    Down hsa-miR-200c-3p 1.7034019   12.256227 9.981328
# 77      5029 -1.547808   4.842417 -8.424155 3.494616e-11 7.926514e-10 15.296172    PRAD  P2RY2    Down  hsa-miR-339-5p 0.6789911    3.208544 5.192101
#     miR_P.Value miR_adj.P.Val    miR_B MirStatus   Estimate   StdError    tValue       pValue     uniHR uniHR_Lower uniHR_Upper unicoxPvalue logrankPvalue
# 64 7.186197e-14  1.882784e-12 21.30487        Up -0.1941886 0.06279724 -3.092311 0.0021144460 0.8238525   0.7075650   0.9592518  0.012566450   0.161731427
# 77 3.224027e-06  1.141480e-05  3.88833        Up -0.4026902 0.11446247 -3.518098 0.0004784352 0.7935574   0.6757902   0.9318475  0.004784937   0.002810577
#    numOfPats numOfEvents unicoxQvalue logrankQvalue             mulHR       mulHR_Lower       mulHR_Upper        mulcoxPvalue            covariates mulcoxQvalue
# 64       479          88   0.03676850    0.30420911 0.824073129144404 0.706983842209402 0.960554515723583  0.0133363232799618 age+histological_type   0.04092737
# 77       479          88   0.02158507    0.01846121  0.79412209893402 0.676038068310838 0.932831947749744 0.00500880171294832 age+histological_type   0.02327620
#      oneYAUC     oneYAUC95CI threeYAUC   threeYAUC95CI CIndex           95%CI
# 64 0.6681761 0.67(0.57-0.76) 0.5772405 0.58(0.49-0.66)    0.6 0.60(0.54-0.66)
# 77 0.6093016 0.61(0.52-0.70) 0.6065000 0.61(0.52-0.69)    0.6 0.60(0.54-0.66)


# nRDiffExp %>% summarise(n_distinct(Symbol), n_distinct(miRNA), n_distinct(DISEASE))
#   n_distinct(Symbol) n_distinct(miRNA) n_distinct(DISEASE)
# 1                 31                51                  14
# nRDiffExp %>% group_by(Symbol, miRNA) %>% distinct(DISEASE, .keep_all = T) %>% summarise(n = n_distinct(DISEASE)) %>% arrange(desc(n)) %>% nrow()
# [1] 61


# nRDiffExp %>% group_by(Symbol, miRNA) %>% distinct(DISEASE, .keep_all = T) %>% summarise(n = n_distinct(DISEASE)) %>%
#   arrange(desc(n)) %>% ungroup() %>% slice_head(n = 5)
#   Symbol miRNA               n
# 1 ADRA2A hsa-miR-23a-3p      3
# 2 ADRA2C hsa-miR-96-5p       3
# 3 ADRA1B hsa-miR-17-5p       2
# 4 CHRNB2 hsa-miR-148b-3p     2
# 5 DRD1   hsa-miR-106b-5p     2


write.table(nRDiffExp, file = 'D:/CancerNeuroscience/Github/result/section4/miRnaNrExprExplain.txt', 
            col.names = T, row.names = F, sep = '\t', quote = FALSE)

########################
mirExpAsso <- nRDiffExp %>% mutate(geneDisease = paste(Symbol, DISEASE, sep = '_'), miRNA = gsub('hsa-', '', miRNA))

geneStatus <- reshape2::dcast(mirExpAsso, miRNA ~ geneDisease, value.var = "NR_logFC") %>% column_to_rownames(var = 'miRNA')
mirStatus <- reshape2::dcast(mirExpAsso, miRNA ~ geneDisease, value.var = "miR_logFC") %>% column_to_rownames(var = 'miRNA')
pfiStatus <- dcast(mirExpAsso, miRNA ~ geneDisease, value.var = "mulHR") %>% column_to_rownames(var = 'miRNA')


gcolFun  <- colorRamp2(breaks = c(-6.0, 0, 5.5), colors = c("#3878C1", "white", "#AB221F"))
mcolFun  <- colorRamp2(breaks = c(-2.0, 0, 4.5), colors = c( "#9FBA95", "white","#B696B6"))


DiagFunc <- function(gStatus, mStatus, surStatus){
  function(j, i, x, y, width, height, fill){
    
    if(!is.na(mStatus[i, j]))
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = gpar(fill = mcolFun(mStatus[i, j]), col = "black"))
    
    if(!is.na(gStatus[i, j]))
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = gcolFun(gStatus[i, j]), col = "black"))    
    
    
    
    if(is.na(surStatus[i, j])){
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(col = "white", fill = NA))
      
      
    }else if(surStatus[i, j] > 1){
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(col = '#FB8C62', fill = NA))
      
    }else{
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(col = "#8C9FCA", fill = NA))  
      
    }
  }
}


mirExpAssoPlot <- Heatmap(geneStatus, rect_gp = gpar(type = "none"), row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                 show_heatmap_legend = F, cluster_rows = F, cluster_columns = F, column_names_rot = 60,
                 cell_fun = DiagFunc(gStatus = geneStatus, mStatus = mirStatus, surStatus = pfiStatus))

lgd <- list(Legend(title = "Exp log2Fc", col_fun = gcolFun, at = c(-5.0, -3.0, 0, 3.0, 5.0), direction = "horizontal"), 
            Legend(title = "miRNA log2Fc", col_fun = mcolFun, at = c(-1.5, -0.5, 0, 2.0, 4.0), direction = "horizontal"))


pdf(file = 'D:/CancerNeuroscience/Github/result/section4/miRNANRExprAssoPlot.pdf', height = 16/2.54, width = 20/2.54)

draw(mirExpAssoPlot, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",heatmap_legend_side = "bottom",merge_legend = TRUE)

dev.off()

