#####################
load(file = '/data/panCanGeneExpData.RData')
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


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

tcgaExpData <- panCanTurGeneExp


nReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
difExpQvalue <- readRDS(file = '/data/lihcCluster1DiffExpGene.rds')
difExpQvalue <- difExpQvalue[difExpQvalue < 0.25]

topGenes <- nReceptors %>% subset.data.frame(Approved.symbol %in% names(difExpQvalue)) %>% remove_rownames() %>% 
  column_to_rownames(var = 'NCBI.Gene.ID') 

disZscoreExp <- tcgaExpData[rownames(topGenes), ]
rownames(disZscoreExp) <- topGenes$Approved.symbol


#####################
hcPearsonWardAverCluster <- readRDS(file = '/result/section5/unFilter/hcPearsonWardAverCluster.rds')
lihcCluster <- hcPearsonWardAverCluster[[5]]$consensusClass %>% as.data.frame() %>% rownames_to_column(var='SAMPLE_BARCODE')
colnames(lihcCluster)[2] <- 'Clusters'

lihcCluster <- lihcCluster %>% left_join(tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE') %>% 
  mutate(Clusters = paste0('Cluster', Clusters)) %>% column_to_rownames(var = 'SAMPLE_BARCODE') %>% subset(DISEASE == 'LIHC')


#####################
# data prepare
disZscoreExp <- disZscoreExp[, intersect(rownames(lihcCluster), colnames(disZscoreExp))] %>% t()
disZscoreExp[is.na(disZscoreExp)] <- 0
lihcCluster <- lihcCluster %>% select(Clusters) %>% mutate(Clusters = ifelse(Clusters == 'Cluster1', 1, 0)) %>% as.matrix()


# Does k-fold cross-validation for glmnet,
set.seed(123)
cv.fit = glmnet::cv.glmnet(x = disZscoreExp, y = lihcCluster, type.measure = 'auc', nfolds=10, family = "binomial")


# plot the cross-validation curve
pdf('/result/section5/lihclike/cv_curve_auc.pdf')
plot(cv.fit)
dev.off()


#####################
# extract non-zero coefficients
est.coef = coef(cv.fit, s = cv.fit$lambda.min)
est.coef <- est.coef[, 's1']
active.k.vals = est.coef[which(est.coef != 0)]
active.k.vals <- active.k.vals[-grep('(Intercept)', names(active.k.vals))]
active.k.vals <- data.frame(symbol = names(active.k.vals), coef = active.k.vals)


est.coef = coef(cv.fit, s = cv.fit$lambda.1se)
est.coef <- est.coef[, 's1']
active.k.vals.1se = est.coef[which(est.coef != 0)]
active.k.vals.1se <- active.k.vals.1se[-grep('(Intercept)', names(active.k.vals.1se))]
active.k.vals.1se <- data.frame(symbol = names(active.k.vals.1se), coef = active.k.vals.1se)


########
geneInfo <- read.csv(file = '/data/gene_with_protein_product.txt', sep = '\t', header = T, stringsAsFactors = F)

active.k.vals <- geneInfo %>% select(symbol, entrez_id, ensembl_gene_id) %>% 
  mutate(entrez_id = as.character(entrez_id)) %>% inner_join(active.k.vals, by = join_by(symbol))

active.k.vals.1se <- geneInfo %>% select(symbol, entrez_id, ensembl_gene_id) %>% 
  mutate(entrez_id = as.character(entrez_id)) %>% inner_join(active.k.vals.1se, by = join_by(symbol))

########

save(cv.fit, active.k.vals, active.k.vals.1se, file='/data/lihc_lasso_binomial_res.RData')


