
library('msigdbr')
library("GSVA")

# load data
tcga.lgg.cli.data <- readRDS(file='D:/CancerNeuroscience/Github/data/LggRiskScores/tcga_lgg_risk_score.rds')
load(file = 'D:/CancerNeuroscience/Github/data/panCanGeneExpData.RData')
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


#####################compare C4-like vs. Others 
# load data
tcga.lgg.cli.data <- readRDS(file='D:/CancerNeuroscience/Github/data/LggRiskScores/tcga_lgg_risk_score.rds')
tcga.lgg.cli.data$risk.categ <- factor(tcga.lgg.cli.data$risk.categ, levels = c('Others', 'C4-like'))

ssgsea_es <- readRDS(file = 'D:/CancerNeuroscience/data/tcga_lgg_hallmark_score.rds')

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

# 
# write.table(p.values, file = 'D:/CancerNeuroscience/Github/result/section5/lgglike/hallmark_score_com_intcga.txt',
#             row.names = T, col.names = T, sep ='\t', quote = F)
# 

##################### plot
hallmarks <- c('HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 
'HALLMARK_P53_PATHWAY', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_COMPLEMENT', 'HALLMARK_INFLAMMATORY_RESPONSE')

plots <- lapply(hallmarks, function(hallmark){
  
  p <- ggboxplot(tcga.lgg.cli.data, x = "risk.categ", y = hallmark, 
                 color = "risk.categ", palette = "jco", add.params = list(size = 0.6), 
                 add = "jitter", xlab = F, ylab = gsub(pattern = 'HALLMARK_', replacement = '', hallmark)) +
    stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", label.x = 1.5) + 
    scale_x_discrete(labels = c('C4-like' = "C4-like", 'Others' = "Others")) + 
    theme(legend.position = "none", axis.title = element_text(size = 8), axis.text = element_text(size = 7),
          legend.title = element_text(size = 5), legend.text = element_text(size = 5))
  
  return(p)
})


ggsave(cowplot::plot_grid(plotlist = plots, ncol=4, nrow=2), width = 10, height = 8, units = 'cm', 
       file='D:/CancerNeuroscience/Github/result/section5/lgglike/select_hallmark_score_com.pdf')

