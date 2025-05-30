

tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')
neurotransmitterReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')

setwd('/data/Taylor_CancerCell_2018_CancerAneuploidy')
copyNumData <- read.csv(file = 'all_thresholded.by_genes_whitelisted.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


copyNumData <- copyNumData %>% column_to_rownames(var = 'Gene.Symbol') %>% select(!c(Locus.ID:Cytoband))
colnames(copyNumData) <- gsub('\\.', '-', substr(colnames(copyNumData), 1, 15))


copyNumData <- copyNumData[, intersect(tcgaPanCanSamples$SAMPLE_BARCODE, colnames(copyNumData))]


copyNumData[copyNumData == -1 | copyNumData == 1] <- 0
nRCopyNumData <- copyNumData[rownames(neurotransmitterReceptors), ]


save(nRCopyNumData, copyNumData, file = '/data/panCanNrCopyNumData.RData')

