
library('msigdbr')
library("GSVA")

# load data
tcga.lgg.cli.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
load(file = '/data/panCanGeneExpData.RData')
tcga.lgg.exp <- panCanTurGeneExp[, intersect(colnames(panCanTurGeneExp), tcga.lgg.cli.data$patient_id)]



# gene set variation analysis
h_gene_sets = msigdbr(species = "human", category = "H")
h_gene_sets = split(x = h_gene_sets$entrez_gene, f = h_gene_sets$gs_name)


## build GSVA parameter object
ssgseapar <- ssgseaParam(
  exprData = as.matrix(tcga.lgg.exp),
  geneSets = h_gene_sets,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE)

## estimate GSVA enrichment scores
ssgsea_es <- gsva(ssgseapar)

# saveRDS(ssgsea_es, file = '/data/tcga_lgg_hallmark_score.rds')


#####################compare high risk vs. low risk 
# load data
tcga.lgg.cli.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
tcga.lgg.cli.data$risk.categ <- factor(tcga.lgg.cli.data$risk.categ, levels = c('low risk', 'high risk'))

ssgsea_es <- readRDS(file = '/data/tcga_lgg_hallmark_score.rds')

tcga.lgg.cli.data <- cbind(tcga.lgg.cli.data, t(ssgsea_es[, tcga.lgg.cli.data$patient_id]))

p.values <- sapply(rownames(ssgsea_es), function(hallmark){
  
  gseascore <- tcga.lgg.cli.data[, c(hallmark, 'risk.categ')]
  colnames(gseascore) <- c('hallmark', 'risk.categ')
  
  meanScore <- gseascore %>% dplyr::select(hallmark, risk.categ) %>% group_by(risk.categ) %>% 
    summarise(meanscore = mean(hallmark))
  
  
  p.value <- wilcox.test(hallmark ~ risk.categ, alternative = "two.sided", data = gseascore)$p.value
  return(setNames(c(meanScore$meanscore, p.value), c('Lsocre', 'Hscore', 'pvalue')))
}, simplify = F)

p.values <- do.call(rbind, p.values) %>% as.data.frame() %>% 
  mutate(qvalue = p.adjust(pvalue, method = 'fdr')) %>% arrange(qvalue) %>% mutate(log2FC = log2(Lsocre/Hscore))


# write.table(p.values, file = '/result/section6/lgglike/hallmark_score_com_intcga.txt',
#             row.names = T, col.names = T, sep ='\t', quote = F)


##################### plot
hallmarks <- c('HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 
'HALLMARK_P53_PATHWAY', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_COMPLEMENT', 'HALLMARK_INFLAMMATORY_RESPONSE')

plots <- lapply(hallmarks, function(hallmark){
  
  p <- ggboxplot(tcga.lgg.cli.data, x = "risk.categ", y = hallmark, 
                 color = "risk.categ", palette = "jco", add.params = list(size = 0.6), 
                 add = "jitter", xlab = F, ylab = hallmark) +
    stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label.x = 1.5) + 
    scale_x_discrete(labels = c('high risk' = "High risk", 'low risk' = "Low risk")) + theme(legend.position = "none")
  
  return(p)
})


ggsave(cowplot::plot_grid(plotlist = plots, ncol=4, nrow=4), 
       file='/result/section6/lgglike/select_hallmark_score_com.pdf')

