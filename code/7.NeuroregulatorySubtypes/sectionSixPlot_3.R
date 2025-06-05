
library('dplyr')
library('tibble')
library('RColorBrewer')
library('pheatmap')
####################

load(file = '/data/panCanGeneExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')


# Z-normalization
panCanTurGeneExp <- as.data.frame(t(panCanTurGeneExp))
panCanTurGeneExp <- panCanTurGeneExp %>% rownames_to_column(var = "SAMPLE_BARCODE")
panCanTurGeneExp <- merge(panCanTurGeneExp, tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE')


panCanTurGeneExp <- split.data.frame(panCanTurGeneExp, f = panCanTurGeneExp$DISEASE)

panCanTurGeneExp <- lapply(panCanTurGeneExp, function(geneExp){
  
  geneExp <- geneExp %>% mutate(DISEASE = NULL) %>% remove_rownames() %>% column_to_rownames(var = 'SAMPLE_BARCODE')
  
  geneExp <- t(scale(geneExp, center = TRUE, scale = TRUE))
  return(geneExp)
})

panCanTurGeneExp <- do.call(cbind, panCanTurGeneExp)

nrPanCanTurGeneExp <- panCanTurGeneExp[as.character(neurotransmitterReceptors$NCBI.Gene.ID), ]
rownames(nrPanCanTurGeneExp) <- neurotransmitterReceptors$Approved.symbol

####################
setwd('/result/section6/consensusCluster_New/New3/')

panCanCluster <- read.csv(file='panCanClusterNew3.k=5.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% rename(SAMPLE_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE') %>% arrange(Clusters, DISEASE) %>% 
  column_to_rownames(var = 'SAMPLE_BARCODE') %>%  mutate(Clusters = paste0('Cluster', Clusters))


rowAnno <- neurotransmitterReceptors %>% select(classesOfNeuroreceptors) %>% rename(nrGroup = classesOfNeuroreceptors)

mycol <- brewer.pal(5, "Set3")
annColors <- list(Clusters = c("Cluster1" = mycol[1], "Cluster2" = mycol[2],
                               "Cluster3" = mycol[3], "Cluster4" = mycol[4], "Cluster5" = mycol[5]))


tmp <- nrPanCanTurGeneExp[, rownames(panCanCluster)]
tmp <- tmp[apply(tmp, 1, function(x) sum(is.na(x))) < 1500, ]
tmp[tmp > 2] <- 2
tmp[tmp < -2] <- -2


# lggCluster <- subset(panCanCluster, DISEASE == 'LGG')
# lggExp <- tmp[, rownames(lggCluster)]

difExpQvalue <- readRDS(file = '/data/ccClusterDiffExpGene_1.rds')
tmp <- tmp[intersect(rownames(tmp), names(head(sort(difExpQvalue), 100))), ]

clusterMaxExpGene <- readRDS(file = '/data/clusterMaxExpGene.rds')

tmp <- tmp[intersect(clusterMaxExpGene$Approved.symbol, rownames(tmp)), ]

symbol <- c('ADRA2C', 'DRD4', 'GABBR1', 'GABRD', 'GRIN2C', 'OPRL1', 
            'ADRA2A', 'ADRB2', 'CHRNA1', 'CHRNA6', 'HRH1', 'HTR2A', 
            'CHRNB2', 'GABBR2', 'GABRB3', 'GLRB', 'GRIK2', 'GRM1', 
            'OPRK1', 'CHRM3', 'GABRA3', 'HTR3C', 'HTR2C', 'GABRR1')

pdf(file = '/result/section6/panCanClusterExp_1.pdf', height = 10)
pheatmap(mat = tmp[symbol, ], cluster_cols = F, cluster_rows = F, treeheight_row = 10, fontsize_row = 6, fontsize = 8, 
         color = colorRampPalette((c("#3878C1", "white","#AB221F")))(100), border_color = NA,
         annotation_col = panCanCluster, annotation_colors = annColors, annotation_row = rowAnno,
         show_colnames = F, show_rownames = T, na_col = "white")

dev.off()




####################
setwd('/result/section6/consensusCluster_New/New3/')

panCanCluster <- read.csv(file='panCanClusterNew3.k=5.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% rename(SAMPLE_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE', 'SUBTYPE')], by = 'SAMPLE_BARCODE') %>% arrange(Clusters, DISEASE) %>% 
  column_to_rownames(var = 'SAMPLE_BARCODE') %>%  mutate(Clusters = paste0('Cluster', Clusters))


panCanCluster <- subset(panCanCluster, DISEASE %in% c('BRCA', 'CESC', 'COAD', 'ESCA', 'GBM', 'HNSC', 'LGG', 'READ', 'SARC', 'STAD', 'TGCT', 'UCEC'))
panCanCluster <- subset(panCanCluster, SUBTYPE != 'NA')

#######
fisherF <- function(data) {
  
  print(table(data[, c('Clusters', 'SUBTYPE')]))
  res <- chisq.test(table(data[, c('Clusters', 'SUBTYPE')]))
  
  return(data.frame(pvalue = res$p.value))

}

enrichPvalue <- panCanCluster %>% group_by(DISEASE) %>% do(fisherF(.))
enrichPvalue$fdr <- p.adjust(enrichPvalue$pvalue, method = 'fdr')

# DISEASE   pvalue      fdr
# <chr>      <dbl>    <dbl>
#   1 BRCA    8.72e-46 4.36e-45
# 2 CESC    7.52e- 3 9.40e- 3
# 3 COAD    1.98e- 1 2.20e- 1
# 4 GBM     5.70e- 3 8.15e- 3
# 5 HNSC    1.02e-18 3.40e-18
# 6 LGG     1.73e-70 1.73e-69
# 7 READ    3.07e- 1 3.07e- 1
# 8 SARC    1.38e- 4 2.75e- 4
# 9 TGCT    8.95e-17 2.24e-16
# 10 UCEC    4.39e- 3 7.32e- 3


# DISEASE   pvalue      fdr
# <chr>      <dbl>    <dbl>
#   1 BRCA    1.91e-96 2.29e-95
# 2 CESC    1.20e- 4 1.81e- 4
# 3 COAD    2.18e- 3 2.38e- 3
# 4 ESCA    1.43e- 8 3.43e- 8
# 5 GBM     1.01e- 3 1.34e- 3
# 6 HNSC    1.41e- 3 1.69e- 3
# 7 LGG     5.21e-72 3.12e-71
# 8 READ    1.54e- 1 1.54e- 1
# 9 SARC    2.27e- 7 4.53e- 7
# 10 STAD    1.36e- 6 2.32e- 6
# 11 TGCT    3.12e-16 1.13e-15
# 12 UCEC    3.78e-16 1.13e-15
#######
panCanCluster <- panCanCluster %>% mutate(Clusters = gsub('Cluster', '', Clusters))


pdf(file = '/result/section6/ccClusterAssoWithMolSubtype_1.pdf')
sapply(split.data.frame(panCanCluster, f = ~DISEASE), function(data){
  

  data <- data %>% arrange(desc(SUBTYPE))
  data$Clusters <- as.numeric(data$Clusters)
  
  pheatmap(mat = data[, 'Clusters', F], cluster_cols = F, cluster_rows = F, treeheight_row = 10, fontsize_row = 6, fontsize = 8, 
           color = brewer.pal(5, "Set3"), border_color = NA,
           annotation_row = data[, 'SUBTYPE', F], show_colnames = F, show_rownames = F, na_col = "black", main = data$DISEASE[1])
})

dev.off()





