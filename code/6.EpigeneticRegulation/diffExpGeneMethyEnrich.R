library('dplyr')
library('tibble')
library('ggplot2')
########################
diffMethy <- readRDS(file = '/data/panCanDiffMethy.rds')


# Diff Expr
nRDiffExp <- readRDS(file = '/data/panCanNrDiffExp.rds')
nRDiffExp <- nRDiffExp %>% mutate(Estatus = ifelse(logFC > 1 & adj.P.Val < 0.05, 'Up', ifelse(logFC < -1 & adj.P.Val < 0.05, 'Down', 'Neutral')))
colnames(nRDiffExp)[9] <- 'Symbol'

nRDiffExp <- nRDiffExp %>% inner_join(diffMethy, by = join_by(DISEASE, Symbol))


# nRDiffExp <- nRDiffExp %>% select(DISEASE, Symbol, Estatus, Mstatus) %>%
#   mutate(Mstatus = factor(Mstatus, levels = c('Hyper', 'Hypo', 'Neutral'))) %>% arrange(Mstatus) %>% distinct(DISEASE, Symbol, Estatus, .keep_all = T)


# > table(nRDiffExp$Estatus, nRDiffExp$Mstatus)
#         Hyper Hypo Neutral
# Down       54    3     213
# Neutral   109    9     595
# Up         17    4     155

# Down
# > fisher.test(matrix(c(54, 216, 109, 604), 2, byrow = T))$p.value
# [1] 0.08378655
# > fisher.test(matrix(c(54, 216, 109, 604), 2, byrow = T))$estimate
# odds ratio
# 1.384852

# Up
# > fisher.test(matrix(c(4, 172, 9, 704), 2, byrow = T))$p.value
# [1] 0.3009801
# > fisher.test(matrix(c(4, 172, 9, 704), 2, byrow = T))$estimate
# odds ratio 
# 1.817674 


########################
dataStat <- nRDiffExp %>% group_by(Estatus, Mstatus) %>% summarise(n = n()) %>% group_by(Estatus) %>% mutate(perc = n/sum(n)) %>% 
  mutate(Mstatus = factor(Mstatus, levels = c('Hyper', 'Neutral', 'Hypo')))

nrDownEnrichHyperPlot <- ggplot(dataStat, aes(x = Estatus, y = perc, fill = Mstatus)) +
  geom_bar(stat = "identity", position = "stack") + labs(y = 'Percentage') + theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5), 
        legend.position = "bottom", legend.justification = "left", axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6)) 


ggsave(filename = /result/section4/nrDownEnrichHyperPlot.pdf', 
       plot = nrDownEnrichHyperPlot, width = 7, height = 7, units = 'cm')

