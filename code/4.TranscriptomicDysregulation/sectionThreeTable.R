

############TCGA
nRDiffExp <- readRDS(file = '/data/panCanNrDiffExp.rds')

nRDiffExp <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% 
  mutate(status = ifelse(logFC > 0,  'Up', 'Down'))


write.table(nRDiffExp, file = '/result/section3/nrDiffExp.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

############Validation
nRDiffExp <- readRDS(file = '/data/curatedPanCanNrDiffExp.rds')

nRDiffExp <- nRDiffExp %>% subset(abs(logFC) > 1 & adj.P.Val < 0.05) %>% 
  mutate(status = ifelse(logFC > 0,  'Up', 'Down'))


write.table(nRDiffExp, file = '/result/section3/ValidationnrDiffExp.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

