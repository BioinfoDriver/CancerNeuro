#####################
load(file = '/data/panCanGeneExpData.RData')
tcgaExpData <- panCanTurGeneExp
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')


neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
difExpQvalue <- readRDS(file = '/data/lggCluster2DiffExpGene.rds')
difExpQvalue <- difExpQvalue[difExpQvalue < 0.05]

topGenes <- neurotransmitterReceptors %>% subset.data.frame(Approved.symbol %in% names(difExpQvalue)) %>% remove_rownames() %>% 
  column_to_rownames(var = 'NCBI.Gene.ID') 

disZscoreExp <- tcgaExpData[rownames(topGenes), ]
rownames(disZscoreExp) <- topGenes$Approved.symbol


#####################
setwd('/result/section6/consensusCluster_New/New3')
panCanCluster <- read.csv(file='panCanClusterNew3.k=5.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
lggCluster <- panCanCluster %>% dplyr::rename(SAMPLE_BARCODE = V1, Clusters = V2) %>% 
  left_join(tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE') %>% 
  arrange(Clusters, DISEASE) %>% mutate(Clusters = paste0('Cluster', Clusters)) %>% 
  column_to_rownames(var = 'SAMPLE_BARCODE') %>% subset(DISEASE == 'LGG')


#####################
# data prepare
disZscoreExp <- disZscoreExp[, intersect(rownames(lggCluster), colnames(disZscoreExp))] %>% t()
disZscoreExp[is.na(disZscoreExp)] <- 0
lggCluster <- lggCluster %>% select(Clusters) %>% mutate(Clusters = ifelse(Clusters == 'Cluster2', 1, 0)) %>% as.matrix()


# Does k-fold cross-validation for glmnet,
set.seed(123)
cv.fit = glmnet::cv.glmnet(x = disZscoreExp, y = lggCluster, type.measure = 'auc', nfolds=10, family = "binomial")


# plot the cross-validation curve
pdf('/result/section6/lgglike/cv_curve_auc.pdf')
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

save(cv.fit, active.k.vals, active.k.vals.1se, file='/data/lgg_lasso_binomial_res.RData')

