
setwd('/data/GDCReferenceFiles')

geneInfo <- read.csv(file = 'gencode.gene.info.v22.tsv', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
geneInfo <- subset(geneInfo, gene_type == 'protein_coding')

mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

geneInfo <- subset(geneInfo, gene_name %in% mc3MutData$Hugo_Symbol)


geneInfo <- geneInfo %>% mutate(breaks = cut(full_length, 
                                             breaks = round(quantile(full_length, probs = seq(0, 1, 0.025))), include.lowest=TRUE))

geneInfo <- geneInfo %>% mutate(exon_breaks = cut(exon_length, 
                                             breaks = round(quantile(exon_length, probs = seq(0, 1, 0.025))), include.lowest=TRUE))


saveRDS(geneInfo, file = '/data/geneInfo.rds')
