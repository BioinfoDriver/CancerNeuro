
library('ConsensusClusterPlus')
library('dplyr')
library('RColorBrewer')
library('pheatmap')


ccRes <- readRDS(file = '/result/section6/consensusCluster_New/New3/ccRes.rds')

##############################
cdfAuc <- sapply(2:9, function(index){
  
  cM <- ccRes[[index]]$consensusMatrix
  Fn = ecdf(cM[lower.tri(cM)])
  
  auc <- integrate(Fn, lower = 0.05, upper = 0.95, subdivisions = 1000L)$value
  return(auc)
})


#########
pdf(file = '/result/section6/deltaArea_1.pdf')

plot(x = 2:9, y = c(cdfAuc[1], diff(cdfAuc)), xlab = 'k', 
     ylab = 'relative change in area under CDF curve', type = 'b', main = 'Delta area')

dev.off()


#########
cols = c('#551F33', '#CBBBC1', '#E4B7BC',  '#BD4146','#ECC68C', '#F5E4CB', '#D6CDBE', '#DDDEDE')


cM <- ccRes[[2]]$consensusMatrix
Fn = ecdf(cM[lower.tri(cM)])

pdf(file = '/result/section6/consensusCDF_1.pdf')

plot(Fn, xlab = 'consensus index', ylab = 'CDF', main = 'consensus CDF', col = cols[1])

for(index in 3:9){
  
  cM <- ccRes[[index]]$consensusMatrix
  Fn = ecdf(cM[lower.tri(cM)])
  
  lines(Fn, col = cols[index-1])

}

legend("bottomright", legend = c(2:9), col = cols, lty = 1)
dev.off()
#########


annCol <- as.data.frame(ccRes[[5]]$consensusClass)
colnames(annCol) <- 'Clusters'
annCol <- annCol %>% mutate(Clusters = paste0('Cluster', Clusters))
  

mycol <- brewer.pal(5, "Set3")
annColors <- list(Clusters = c("Cluster1" = mycol[1], "Cluster2" = mycol[2],
                              "Cluster3" = mycol[3], "Cluster4" = mycol[4], "Cluster5" = mycol[5]))

heatdata <- ccRes[[5]]$consensusMatrix
dimnames(heatdata) <- list(rownames(annCol), rownames(annCol))



ccplot <- pheatmap(mat = heatdata, cluster_cols = T, cluster_rows = T,
         color = colorRampPalette((c("white","steelblue")))(100), border_color = NA,
         annotation_col = annCol, annotation_colors = annColors,
         show_colnames = F, show_rownames = F)

pdf(file = '/result/section6/consensusMatrixCluster5_1.pdf')

dev.off()
#########

