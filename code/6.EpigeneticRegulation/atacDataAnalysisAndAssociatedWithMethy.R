library('GenomicRanges')
library('ggplot2')
library('rstatix')
library('ggpubr')

# ACAT data
tcgaPanCanSamples <- readRDS(file = 'D:/CancerNeuroscience/Github/data/tcgaPanCanSamples.rds')
atacPeakData <- read.csv(file = 'F:/PostdoctoralDataBackup/DesktopCP/F/BigData/UCSCXena/ATAT-seq_Hub/Pan-Cancer(PANCAN)/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz', 
                         header = T, sep = '\t')
colnames(atacPeakData) <- gsub('\\.', '-', colnames(atacPeakData))
colnames(atacPeakData) <- stringr::str_sub(colnames(atacPeakData), 1, 15)

# subset(tcgaPanCanSamples, SAMPLE_BARCODE %in% colnames(atacPeakData)[-1]) %>% group_by(DISEASE) %>% summarise(patsNum = n()) %>% as.data.frame()
# table(substr(colnames(atacPeakData)[-1], 14, 15))
#  01  02  05  06 
# 390   3   1  10 

atacPeakData <- tidyr::pivot_longer(data = atacPeakData, names_to = "SAMPLE_BARCODE", values_to = "ATAC_Score", cols = starts_with("TCGA"))
colnames(atacPeakData)[1] <- 'id'

tcgaPanCanSamples <- readRDS(file = 'D:/CancerNeuroscience/Github/data/tcgaPanCanSamples.rds')

atacPeakData <- inner_join(atacPeakData, tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'DISEASE')], by = 'SAMPLE_BARCODE')


# ACAT peak annotation
atacPeakAnno <- read.csv(file = 'F:/PostdoctoralDataBackup/DesktopCP/F/BigData/UCSCXena/ATAT-seq_Hub/Pan-Cancer(PANCAN)/TCGA_ATAC_peak.all.probeMap', 
         header = T, sep = '\t')
diffMethy <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanDiffMethy.rds')


# Finding overlapping genomic ranges
queryRegions <- diffMethy %>% mutate(end = pos) %>% 
  makeGRangesFromDataFrame(seqnames.field = 'chr', start.field="pos", 
                           end.field='end', strand.field="strand", keep.extra.columns=TRUE)

subjectRegions <- atacPeakAnno %>% makeGRangesFromDataFrame(seqnames.field = 'chrom', start.field="chromStart", 
                                                          end.field='chromEnd', strand.field="strand", keep.extra.columns=TRUE)
subjectRegions <- resize(subjectRegions, width = width(subjectRegions) + 500, fix = "center")


regionOverlaps  <- findOverlaps(queryRegions, subjectRegions, type='within', select = 'all', ignore.strand=TRUE)

qRegions <- queryRegions[queryHits(regionOverlaps)]
sRegions <- subjectRegions[subjectHits(regionOverlaps), ]

qRegions <- as.data.frame(qRegions) %>% dplyr::rename(seqnames_M = seqnames, start_M = start, end_M = end, width_M = width, strand_M = strand)
sRegions <- as.data.frame(sRegions) %>% dplyr::rename(seqnames_A = seqnames, start_A = start, end_A = end, width_A = width, strand_A = strand)

ovRegions <- cbind.data.frame(qRegions, sRegions)


# Finding overlapping genomic ranges
mProbeAtacPeakData <- inner_join(atacPeakData, ovRegions, by = c('id', 'DISEASE'))

mProbeAtacPeakData <- mProbeAtacPeakData %>% mutate(Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Hypo')))
statTest <- mProbeAtacPeakData %>% wilcox_test(ATAC_Score ~ Mstatus) %>% 
  adjust_pvalue(method = "BH") %>% add_xy_position(x = "Mstatus")
# 1.2405e-12 2.8110e-14 1.1400e-09

mProbeAtacPlot <- ggplot(mProbeAtacPeakData, aes(x = Mstatus, y = ATAC_Score)) +
  geom_boxplot(outlier.size = 0.5) + geom_violin(aes(fill = Mstatus), alpha = 0.8) + 
  stat_pvalue_manual(statTest, label = "p.adj", tip.length = 0.01, bracket.size = 0.6, label.size = 3) + 
  scale_fill_manual(values = c('Hyper' = "#F89999", 'Hypo' = "#C77CFF", 'Neutral' = '#99ADFF'), name = NULL) +  
  labs(x = NULL, y = "ATAC z-score") + theme_bw() + 
  theme(legend.position = "top", legend.justification = "left", axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))


# length(unique(mProbeAtacPeakData$Name))
# [1] 31
# length(unique(diffMethy$Name))
# [1] 157

######### Different expression
nRDiffExp <- readRDS(file = 'D:/CancerNeuroscience/Github/data/panCanNrDiffExp.rds')
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(logFC > 1 & adj.P.Val < 0.05, 'Up', 
                                                   ifelse(logFC < -1 & adj.P.Val < 0.05, 'Down', 'Neutral'))) %>% 
  dplyr::rename(Symbol = Approved.symbol) %>% select(DISEASE, Symbol, Estatus)

mProbeAtacPeakData <- mProbeAtacPeakData %>% left_join(nRDiffExp, by = c('Symbol', 'DISEASE'))


######## Plot
nRDiffExpAtacData <- mProbeAtacPeakData %>% subset(!is.na(Estatus))

statTest <- nRDiffExpAtacData %>% wilcox_test(ATAC_Score ~ Estatus) %>% 
  adjust_pvalue(method = "BH") %>% add_xy_position(x = "Estatus")
# 1.824e-46 2.400e-01 1.071e-31

nRDiffExpAtacPlot <- ggplot(nRDiffExpAtacData, aes(x = Estatus, y = ATAC_Score)) +
  geom_boxplot(outlier.size = 0.5) + geom_violin(aes(fill = Estatus), alpha = 0.8) + 
  stat_pvalue_manual(statTest, label = "p.adj", tip.length = 0.01, bracket.size = 0.6, label.size = 3) + 
  scale_fill_manual(values = c('Down' = "#F89999", 'Up' = "#C77CFF", 'Neutral' = '#99ADFF'), name = NULL) +  
  labs(x = NULL, y = "ATAC z-score") + theme_bw() + 
  theme(legend.position = "top", legend.justification = "left", axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))


nRDiffExpFilterAtacData <- mProbeAtacPeakData %>% subset((Estatus == 'Down' & Mstatus == 'Hyper') | (Estatus == 'Up' & Mstatus == 'Hypo') |
                                                     (Estatus == 'Neutral' & Mstatus == 'Neutral'))

statTest <- nRDiffExpFilterAtacData %>% wilcox_test(ATAC_Score ~ Estatus) %>% 
  adjust_pvalue(method = "BH") %>% add_xy_position(x = "Estatus")
# 2.301e-22 4.890e-13 5.900e-06

nRDiffExpFilterAtacPlot <- ggplot(nRDiffExpFilterAtacData, aes(x = Estatus, y = ATAC_Score)) +
  geom_boxplot(outlier.size = 0.5) + geom_violin(aes(fill = Estatus), alpha = 0.8) + 
  stat_pvalue_manual(statTest, label = "p.adj", tip.length = 0.01, bracket.size = 0.6, label.size = 3) + 
  scale_fill_manual(values = c('Down' = "#F89999", 'Up' = "#C77CFF", 'Neutral' = '#99ADFF'), name = NULL) +  
  labs(x = NULL, y = "ATAC z-score") + theme_bw() + 
  theme(legend.position = "top", legend.justification = "left", axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))


ggsave(plot = mProbeAtacPlot + nRDiffExpAtacPlot + nRDiffExpFilterAtacPlot, width = 20, height = 7, units = 'cm',
       filename = 'D:/CancerNeuroscience/Github/result/section4/mProbeAtacScorePlot.pdf')


######################### ATAC vs. Methylation vs. Expression
load(file = 'D:/CancerNeuroscience/Github/data/panCanMethyData.RData')
# panCanTurMethy, panCanPairdTurMethy, panCanPairdNormMethy
load(file = 'D:/CancerNeuroscience/Github/data/panCanGeneExpData.RData')
nReceptors <- readRDS(file = 'D:/CancerNeuroscience/Github/data/neurotransmitterReceptors.rds')
nReceptors <- nReceptors %>% mutate(NCBI.Gene.ID = as.character(NCBI.Gene.ID))

# Expression
nrExp <- panCanTurGeneExp %>% select(any_of(unique(mProbeAtacPeakData$SAMPLE_BARCODE)))
nrExp <- nrExp[nReceptors$NCBI.Gene.ID, ]
rownames(nrExp) <- nReceptors$Approved.symbol
nrExp <- nrExp %>% rownames_to_column(var = 'Symbol')
nrExp <- tidyr::pivot_longer(data = nrExp, names_to = "SAMPLE_BARCODE", values_to = "geneExp", cols = starts_with("TCGA"))

# Methylation
nrMethy <- panCanTurMethy %>% select(any_of(unique(mProbeAtacPeakData$SAMPLE_BARCODE)))
nrMethy <- nrMethy[unique(mProbeAtacPeakData$Name), ]
nrMethy <- nrMethy %>% rownames_to_column(var = 'Name')
nrMethy <- tidyr::pivot_longer(data = nrMethy, names_to = "SAMPLE_BARCODE", values_to = "geneMethy", cols = starts_with("TCGA"))


mProbeAtacPeakData <- mProbeAtacPeakData %>% inner_join(nrMethy, by = c('Name', 'SAMPLE_BARCODE'))
mProbeAtacPeakData <- mProbeAtacPeakData %>% inner_join(nrExp, by = c('Symbol', 'SAMPLE_BARCODE'))


methyExpAcat <- mProbeAtacPeakData %>% select(Name, Symbol, DISEASE, id, geneExp, geneMethy, ATAC_Score, Mstatus, Estatus) %>% 
  distinct() %>% split.data.frame(f = ~ Name + Symbol + DISEASE + id)

methyExpAcat <- methyExpAcat[sapply(methyExpAcat, nrow) > 0]


methyExpAcatCor <- lapply(seq(length(methyExpAcat)), function(i){

  print(i)
  assoData <- methyExpAcat[[i]]
  
  if(sum(is.na(assoData$geneExp))/nrow(assoData) > 0.2){
    res <- data.frame(Name = assoData$Name[1], Symbol = assoData$Symbol[1], DISEASE = assoData$DISEASE[1], id = assoData$id[1], nSams = nrow(assoData), 
                      Mstatus = assoData$Mstatus[1], Estatus = assoData$Estatus[1], emCor = NA, emPvalue = NA, eaCor = NA, eaPvalue = NA, maCor = NA, maPvalue = NA) 
    
  }else{
    emCor <- cor.test(formula = ~ geneExp + geneMethy, data = assoData, method = "spearman", alternative = 'less')$estimate
    emPvalue <- cor.test(formula = ~ geneExp + geneMethy, data = assoData, method = "spearman", alternative = 'less')$p.value
    
    eaCor <- cor.test(formula = ~ geneExp + ATAC_Score, data = assoData, method = "spearman", alternative = 'greater')$estimate
    eaPvalue <- cor.test(formula = ~ geneExp + ATAC_Score, data = assoData, method = "spearman", alternative = 'greater')$p.value  
    
    maCor <- cor.test(formula = ~ geneMethy + ATAC_Score, data = assoData, method = "spearman", alternative = 'less')$estimate
    maPvalue <- cor.test(formula = ~ geneMethy + ATAC_Score, data = assoData, method = "spearman", alternative = 'less')$p.value
    
    
    res <- data.frame(Name = assoData$Name[1], Symbol = assoData$Symbol[1], DISEASE = assoData$DISEASE[1], id = assoData$id[1], nSams = nrow(assoData),
                      Mstatus = assoData$Mstatus[1], Estatus = assoData$Estatus[1], emCor = emCor, emPvalue = emPvalue, eaCor = eaCor, eaPvalue = eaPvalue, maCor = maCor, maPvalue = maPvalue)    
  }
  return(res)
})


methyExpAcatCor <- do.call(rbind, methyExpAcatCor)

methyExpAcatCor %>% subset(Mstatus == 'Hyper' & Estatus == 'Down' & emCor < -0.05 & eaCor > 0.05 & maCor < -0.05)


