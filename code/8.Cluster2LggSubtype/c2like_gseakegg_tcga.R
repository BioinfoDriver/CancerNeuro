
library('ggpubr')
library('ggplot2')
library('limma')
library('clusterProfiler')
library('org.Hs.eg.db')
library('enrichplot')

# load data
tcga.lgg.cli.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
load(file = '/data/panCanGeneExpData.RData')
tcga.lgg.exp <- panCanTurGeneExp[, intersect(colnames(panCanTurGeneExp), tcga.lgg.cli.data$patient_id)]


# diff gene
tcga.lgg.cli.data$risk.categ <- factor(tcga.lgg.cli.data$risk.categ, levels = c('low risk', 'high risk'))
design <- model.matrix(~risk.categ, data = tcga.lgg.cli.data)


fit <- lmFit(tcga.lgg.exp, design)
fit <- eBayes(fit)
difRes <- topTable(fit, coef="risk.categhigh risk", adjust.method = "BH", n = Inf)



# GO,KEGG pathway over-representation analysis
diff.genes <- rownames(subset(difRes, abs(logFC) > 1.0 & adj.P.Val < 0.05))

ego <- enrichGO(gene = diff.genes, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)

ekegg <- enrichKEGG(gene = diff.genes, keyType = "ncbi-geneid", pAdjustMethod = "BH", organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


# KEGG Gene Set Enrichment Analysis
gene.list <- setNames(difRes$logFC, rownames(difRes))
gene.list = sort(gene.list, decreasing = TRUE)

gseakegg <- gseKEGG(geneList = gene.list, organism = 'hsa',
               keyType = "ncbi-geneid", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
               pAdjustMethod = "BH", verbose = FALSE)

gseakegg <- setReadable(gseakegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# save(difRes, ego, ekegg, gseakegg, file = '/data/tcga_lgg_c2like_gokegg_enrich.RData')


################ plot
# write.table(ego, file = '/result/section6/lgglike/goenrich.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# write.table(ekegg, file = '/result/section6/lgglike/keggenrich.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# write.table(gseakegg, file = '/result/section6/lgglike/kegggsea.txt', sep = '\t', col.names = T, row.names = F, quote = F)


categorys <- c('Glutamatergic synapse', 'GABAergic synapse', 'Dopaminergic synapse', 'Serotonergic synapse', 'Cholinergic synapse', 
               'GnRH secretion', 
               'Synaptic vesicle cycle', 'Neuroactive ligand-receptor interaction', 
               'cAMP signaling pathway', 'Calcium signaling pathway', 
               'JAK-STAT signaling pathway',
               'NF-kappa B signaling pathway',
               'ECM-receptor interaction',
               'TNF signaling pathway',
               'Focal adhesion',
               'Cell adhesion molecules',
               'Cytokine-cytokine receptor interaction', 'Th17 cell differentiation', 
               'Th1 and Th2 cell differentiation',
               'Antigen processing and presentation',
               'Natural killer cell mediated cytotoxicity',
               'Leukocyte transendothelial migration',
               'PI3K-Akt signaling pathway',
               'p53 signaling pathway',
               'Cell cycle',
               'IL-17 signaling pathway',
               'Chemokine signaling pathway')


pdf('/result/section6/lgglike/gsea_kegg.pdf')
dotplot(gseakegg, showCategory=categorys, split=".sign") + facet_grid(.~.sign)

gseaplot2(gseakegg, geneSetID = c(46, 61, 63), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

gseaplot2(gseakegg, geneSetID = c(6, 7, 41), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

gseaplot2(gseakegg, geneSetID = c(33, 48, 89), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

gseaplot2(gseakegg, geneSetID = c(18, 24, 53), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")

dev.off()
