library('dplyr')
library('tibble')
library('reshape2')
library('ComplexHeatmap')
library('circlize')

#############
diffMethy <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanDiffMethy.rds')
nRDiffExp <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanNrDiffExp.rds')
panCanSurPvalue <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanRhoSurUniMulLogPvalue_OS.rds')

#############
diffMethy <- diffMethy %>% subset(Mstatus != 'Neutral') %>% dplyr::rename(MlogFC = logFC)

nRDiffExp <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val <= 0.05) %>% 
  mutate(Estatus = ifelse(logFC > 0,  'Up', 'Down')) %>% dplyr::rename(Symbol = Approved.symbol)

MEassoRes <- nRDiffExp %>% inner_join(diffMethy, by = c('DISEASE', 'Symbol')) %>% 
  subset((Mstatus == 'Hyper' & Estatus == 'Down') | (Mstatus == 'Hypo' & Estatus == 'Up'))


panCanSurPvalue <- subset(panCanSurPvalue, mulcoxQvalue <= 0.1)
panCanSurPvalue <- panCanSurPvalue %>% mutate(Direction = ifelse(mulHR > 1.0, 'High level-> worse survival', 'High level-> better survival')) %>% 
  dplyr::rename(DISEASE = disease, Symbol = geneSymbol)


MEassoRes <- MEassoRes %>% left_join(panCanSurPvalue, by = join_by(DISEASE, Symbol))


na.omit(MEassoRes) %>% subset((Estatus == 'Down' & Direction == 'High level-> better survival') | 
                                (Mstatus == 'Up' & Estatus == 'High level-> worse survival'))

##################OS
#    geneID     logFC  AveExpr         t     P.Value    adj.P.Val        B DISEASE Symbol Estatus  chr       pos strand       Name Enhancer RefGene_Group   MlogFC
# 24   2890 -5.011651 6.110474 -19.52144 4.08642e-27 1.121441e-24 51.51262    LUAD  GRIA1    Down chr5 152870490      - cg08578734     TRUE       1stExon 1.022354
# 25   2890 -5.011651 6.110474 -19.52144 4.08642e-27 1.121441e-24 51.51262    LUAD  GRIA1    Down chr5 152870258      + cg17020834     TRUE       1stExon 1.409823
# 26   2890 -5.011651 6.110474 -19.52144 4.08642e-27 1.121441e-24 51.51262    LUAD  GRIA1    Down chr5 152870258      + cg17020834     TRUE         5'UTR 1.409823
#        mDiff       pValue      adjPVal Mstatus     uniHR  uniHR_Lower uniHR_Upper unicoxPvalue logrankPvalue numOfPats numOfEvents unicoxQvalue logrankQvalue
# 24 0.2492013 3.632395e-07 1.222440e-06   Hyper 0.8686236 0.0002271546   0.8059589    0.9361606    0.06953032       493         176            1     0.2275538
# 25 0.2161321 9.806183e-09 7.800373e-08   Hyper 0.8686236 0.0002271546   0.8059589    0.9361606    0.06953032       493         176            1     0.2275538
# 26 0.2161321 9.806183e-09 7.800373e-08   Hyper 0.8686236 0.0002271546   0.8059589    0.9361606    0.06953032       493         176            1     0.2275538
#                mulHR       mulHR_Lower       mulHR_Upper        mulcoxPvalue                 covariates mulcoxQvalue   oneYAUC     oneYAUC95CI threeYAUC
# 24 0.889855796149966 0.816592301127468 0.969692387312975 0.00776723966913063 age+gender+race+ajcc_stage   0.07975237 0.6243454 0.62(0.54-0.71) 0.5960563
# 25 0.889855796149966 0.816592301127468 0.969692387312975 0.00776723966913063 age+gender+race+ajcc_stage   0.07975237 0.6243454 0.62(0.54-0.71) 0.5960563
# 26 0.889855796149966 0.816592301127468 0.969692387312975 0.00776723966913063 age+gender+race+ajcc_stage   0.07975237 0.6243454 0.62(0.54-0.71) 0.5960563
#      threeYAUC95CI CIndex           95%CI                    Direction
# 24 0.60(0.53-0.67)   0.59 0.59(0.54-0.64) High level-> better survival
# 25 0.60(0.53-0.67)   0.59 0.59(0.54-0.64) High level-> better survival
# 26 0.60(0.53-0.67)   0.59 0.59(0.54-0.64) High level-> better survival

##################PFI
# geneID     logFC  AveExpr         t      P.Value    adj.P.Val        B DISEASE Symbol Estatus  chr       pos strand       Name Enhancer RefGene_Group   MlogFC
# 6     148 -1.634161 6.035829 -6.889052 8.643006e-09 7.956770e-08 9.865704    PRAD ADRA1A    Down chr8  26723365      + cg22461835                TSS1500 1.119519
# 35   2892 -1.418260 4.356839 -5.815567 4.135851e-07 2.355966e-06 6.068061    PRAD  GRIA3    Down chrX 122318040      - cg19994834                 TSS200 1.682751
# 36   2892 -1.418260 4.356839 -5.815567 4.135851e-07 2.355966e-06 6.068061    PRAD  GRIA3    Down chrX 122318339      + cg23424962                  5'UTR 2.069922
#        mDiff       pValue      adjPVal Mstatus     uniHR uniHR_Lower uniHR_Upper unicoxPvalue logrankPvalue numOfPats numOfEvents unicoxQvalue logrankQvalue
# 6  0.1749240 1.810088e-09 7.919133e-08   Hyper 0.7529494   0.6655009   0.8518890 6.643635e-06  0.0003110807       479          88  0.000174949    0.01846121
# 35 0.2524787 4.921243e-07 2.100531e-06   Hyper 0.7952198   0.6920533   0.9137656 1.229414e-03  0.0051975142       479          88  0.009712375    0.02645765
# 36 0.2023762 9.896734e-08 7.216369e-07   Hyper 0.7952198   0.6920533   0.9137656 1.229414e-03  0.0051975142       479          88  0.009712375    0.02645765
#                mulHR       mulHR_Lower       mulHR_Upper         mulcoxPvalue            covariates mulcoxQvalue   oneYAUC     oneYAUC95CI threeYAUC   threeYAUC95CI
# 6  0.754561291490832 0.666886465512749 0.853762629863161 7.86857137983011e-06 age+histological_type 0.0002072057 0.6539412 0.65(0.56-0.75) 0.6557450 0.66(0.57-0.74)
# 35 0.800631558131319 0.696271429298828 0.920633627781203  0.00180572502917619 age+histological_type 0.0138991034 0.6011283 0.60(0.49-0.71) 0.6052196 0.61(0.52-0.69)
# 36 0.800631558131319 0.696271429298828 0.920633627781203  0.00180572502917619 age+histological_type 0.0138991034 0.6011283 0.60(0.49-0.71) 0.6052196 0.61(0.52-0.69)
#    CIndex           95%CI                    Direction
# 6    0.65 0.65(0.59-0.71) High level-> better survival
# 35   0.60 0.60(0.54-0.67) High level-> better survival
# 36   0.60 0.60(0.54-0.67) High level-> better survival

# MEassoRes %>% summarise(n_distinct(Symbol), n_distinct(Name), n_distinct(DISEASE))
#   n_distinct(Symbol) n_distinct(Name) n_distinct(DISEASE)
# 1                 26               32                  13

MEassoRes %>% group_by(Symbol, Estatus, Mstatus) %>% distinct(DISEASE, .keep_all = T) %>% summarise(n = n_distinct(DISEASE)) %>%
  arrange(desc(n)) %>% ungroup() %>% slice_head(n = 10)
# Symbol Estatus Mstatus     n
# <chr>  <chr>   <chr>   <int>
# 1 GRIA3  Down    Hyper       4
# 2 ADRA1D Down    Hyper       3
# 3 ADRA2B Down    Hyper       3
# 4 GRIA1  Down    Hyper       3
# 5 GRIN2A Down    Hyper       3
# 6 ADRA1A Down    Hyper       2
# 7 DRD1   Down    Hyper       2
# 8 DRD2   Down    Hyper       2
# 9 GRIK2  Down    Hyper       2
# 10 HTR7  Down    Hyper       2

write.table(MEassoRes, file = 'D:/CancerNeuroscience/Github/result/section4/methyNrExprExplain.txt',
            col.names = T, row.names = F, sep = '\t', quote = F)


########################
geneStatus <- reshape2::dcast(MEassoRes, Symbol ~ DISEASE, fun.aggregate = mean, value.var = "logFC") %>% column_to_rownames(var = 'Symbol')
methyStatus <- reshape2::dcast(MEassoRes, Symbol ~ DISEASE, fun.aggregate = mean, value.var = "MlogFC") %>% column_to_rownames(var = 'Symbol')
osStatus <- reshape2::dcast(MEassoRes, Symbol ~ DISEASE, fun.aggregate = mean, value.var = "mulHR") %>% column_to_rownames(var = 'Symbol')

gcolFun  <- colorRamp2(breaks = c(-5.5, 0, 2.5), colors = c("#3878C1", "white", "#AB221F"))
mcolFun  <- colorRamp2(breaks = c(-2, 0, 3.0), colors = c( "#9FBA95", "white","#B696B6"))


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



methyExpAssoPlot <- Heatmap(geneStatus, rect_gp = gpar(type = "none"), row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                          show_heatmap_legend = F, cluster_rows = F, cluster_columns = F, column_names_rot = 60,
                          cell_fun = DiagFunc(gStatus = geneStatus, mStatus = methyStatus, surStatus = osStatus))


lgd <- list(Legend(title = "Exp log2Fc", col_fun = gcolFun, at = c(-5.5, -3.5, 0, 1.5, 2.5), direction = "horizontal"), 
            Legend(title = "Methy log2Fc", col_fun = mcolFun, at = c(-2, -1, 0, 1.5, 3.0), direction = "horizontal"))


pdf(file = 'D:/CancerNeuroscience/Github/result/section4/methyNRExprAssoPlot.pdf', height = 10/2.54, width = 10/2.54)

draw(methyExpAssoPlot, annotation_legend_list = lgd, 
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom", merge_legend = TRUE)

dev.off()

