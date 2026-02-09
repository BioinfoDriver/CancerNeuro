
library(ggplot2)

###############
nROfmiRNAflm <- readRDS(file = 'D:/CancerNeuroscience/Github/data/nROfmiRNARegulatorFLM.rds')

nROfmiRNAflm <- subset(nROfmiRNAflm, pValue < 0.05)

coefDensityPlot <- ggplot(nROfmiRNAflm, aes(x = Estimate)) +geom_density(fill = "skyblue") + 
  labs(x = "Coefficients", y = "Density") + scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2)) + 
  annotate("text", x = -0.5, y = 1.20, label = "-0.15", color = 'blue') + 
  geom_vline(xintercept = -0.15, linetype = 'dashed', color = 'red', size = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))

ggsave(filename = 'D:/CancerNeuroscience/Github/result/section4/coefDensityPlot.pdf', 
       plot = coefDensityPlot, width = 7, height = 7, units = 'cm')

nROfmiRNAflm <- subset(nROfmiRNAflm, pValue < 0.05 & Estimate < -0.15)
########## Stat
# > dim(nROfmiRNAflm)
# [1] 1062    7


nROfmiRNAflm %>% select(miRNA, nRgene) %>% distinct() %>% nrow()
# 435

nROfmiRNAflm %>% select(miRNA, nRgene) %>% distinct() %>% summarise(n_distinct(miRNA), n_distinct(nRgene))
#   n_distinct(miRNA) n_distinct(nRgene)
# 1               155                 65

nROfmiRNAflm %>% select(miRNA, nRgene) %>% distinct(miRNA, nRgene) %>% 
  count(miRNA) %>% arrange(desc(n)) %>% subset(n >= 2) %>% nrow()
# [1] 100

nROfmiRNAflm %>% select(miRNA, nRgene) %>% distinct(miRNA, nRgene) %>% 
  count(miRNA) %>% arrange(desc(n)) %>% slice_head(n = 20)

#              miRNA  n
# 1  hsa-miR-181b-5p 10
# 2   hsa-miR-15a-5p  9
# 3   hsa-miR-30e-5p  9
# 4    hsa-let-7d-5p  8
# 5    hsa-let-7g-5p  8
# 6   hsa-miR-15b-5p  8
# 7    hsa-miR-16-5p  8
# 8  hsa-miR-181a-5p  8
# 9   hsa-miR-27b-3p  8
# 10  hsa-miR-30b-5p  8
# 11   hsa-let-7e-5p  7
# 12  hsa-miR-27a-3p  7
# 13  hsa-miR-30c-5p  7
# 14  hsa-miR-30d-5p  7
# 15  hsa-miR-424-5p  7
# 16  hsa-miR-497-5p  7
# 17   hsa-let-7b-5p  6
# 18   hsa-let-7i-5p  6
# 19 hsa-miR-181d-5p  6
# 20  hsa-miR-92a-3p  6

nROfmiRNAflm %>% select(miRNA, nRgene) %>% distinct(miRNA, nRgene) %>% count(nRgene) %>% 
  arrange(desc(n)) %>% subset(n >= 2) %>% nrow()
# [1] 59

nROfmiRNAflm %>% select(miRNA, nRgene) %>% distinct(miRNA, nRgene) %>% count(nRgene) %>% 
  arrange(desc(n)) %>% slice_head(n = 20)

#    nRgene  n
# 1  GRIN2A 28
# 2   ADRB1 20
# 3   CHRM3 19
# 4  GABBR2 18
# 5   P2RY2 18
# 6   ADRB2 13
# 7  GABRB2 12
# 8   GRIK2 12
# 9   GRIK3 12
# 10 ADRA2A 11
# 11  GRIA1 11
# 12  GRID1 11
# 13 CHRNA7  9
# 14   DRD1  9
# 15 GABRB3  9
# 16  GRIN1  9
# 17  HTR2A  9
# 18  P2RX1  9
# 19 GABBR1  8
# 20  GABRE  8

nROfmiRNAflm %>% select(miRNA, nRgene, Disease) %>% group_by(miRNA, nRgene) %>% distinct(Disease, .keep_all = T) %>% 
  summarise(n = n()) %>% arrange(desc(n)) %>% subset(n >=2) %>% nrow()
# [1] 242

nROfmiRNAflm %>% select(miRNA, nRgene, Disease) %>% group_by(miRNA, nRgene) %>% distinct(Disease, .keep_all = T) %>% 
  summarise(n = n()) %>% arrange(desc(n)) %>% subset(n >=2) %>% ungroup() %>% slice_head(n = 20)
#     miRNA           nRgene      n
#   1 hsa-let-7g-5p   ADRB2      15
# 2 hsa-miR-182-5p  GABBR1     15
# 3 hsa-miR-148b-3p CHRNB2     13
# 4 hsa-miR-15b-5p  GABBR1     12
# 5 hsa-miR-15b-5p  GRIN1      11
# 6 hsa-miR-23a-3p  ADRA2A     11
# 7 hsa-miR-16-5p   GRIN1      10
# 8 hsa-miR-29a-3p  P2RX5      10
# 9 hsa-miR-30b-5p  ADRA2A      9
# 10 hsa-miR-96-5p   ADRA2C      9
# 11 hsa-miR-24-3p   HRH1        8
# 12 hsa-miR-27a-3p  GABRP       8
# 13 hsa-miR-98-5p   P2RX1       8
# 14 hsa-let-7c-5p   ADRB2       7
# 15 hsa-miR-128-3p  ADORA2B     7
# 16 hsa-miR-130b-3p GRIK2       7
# 17 hsa-miR-15b-5p  ADORA2A     7
# 18 hsa-miR-200a-3p P2RY1       7
# 19 hsa-miR-330-5p  GRIN1       7
# 20 hsa-miR-34a-5p  GRID1       7


write.table(nROfmiRNAflm, file = 'D:/CancerNeuroscience/Github/result/section4/miRnaTarByDisease.txt', 
            col.names = T, row.names = F, sep = '\t', quote = FALSE)

