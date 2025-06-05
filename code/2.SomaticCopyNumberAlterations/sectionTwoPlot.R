
library(dplyr)
library(ggplot2)


nrGisticPeakQvalue <- readRDS(file = '/data/nrGisticPeakQvalue.rds')

##############
nrPeakQvalue <- as.data.frame(nrGisticPeakQvalue)
nrPeakQvalue <- subset(nrPeakQvalue, !(Disease %in% c('COADREAD', 'GBMLGG', 'KIPAN', 'STES')))


nrPeakQvalue$sizeLable <- cut(nrPeakQvalue$q.values, breaks = c(0, 1.0e-05, 1.0e-03, 1.0e-01, 0.25), 
                        labels = c('q < 1.0e-05', 'q < 1.0e-03', 'q < 1.0e-01', 'q < 0.25'))

nrPeakStat <- nrPeakQvalue %>% group_by(Approved.symbol, Direction) %>% count(name = 'numOfDis') %>% 
  reshape2::dcast(Approved.symbol~Direction, value.var = "numOfDis", fill = 0) %>% arrange(desc(Amplification), Deletion)

nrPeakStat <- nrPeakStat %>% mutate(numOfDis = Amplification + Deletion)%>% 
  arrange(desc(Amplification), Deletion) %>% slice(c(1:12, (n()-11):n()))


nrPeakQvalue <- nrPeakQvalue  %>% inner_join(nrPeakStat, by = 'Approved.symbol') %>% # arrange(desc(numOfDis)) %>%
  mutate(Approved.symbol = factor(Approved.symbol, levels= unique(Approved.symbol)))


fillLabel <- c("red", "blue")
names(fillLabel) <- c('Amplification', 'Deletion')

sizeLable <- c(4, 3, 2, 1)
names(sizeLable) <- c('q < 1.0e-05', 'q < 1.0e-03', 'q < 1.0e-01', 'q < 0.25')


plot1 <- ggplot(nrPeakQvalue, aes(x = Disease, y = Approved.symbol, color = Direction, size = sizeLable)) +
  geom_point() + labs(x = "Gene", y = "Cancer Types") + 
  scale_size_manual(values = sizeLable, guide = "legend", name = "Gastic q value") + 
  scale_color_manual(values = fillLabel, guide = "legend", name = "SCNA levels") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# plots1 <- ggplot(subset(nrPeakQvalue, numOfDis > 5), 
#                 aes(x = Disease, y = reorder(Approved.symbol, -index), color = Direction, size = sizeLable)) +
#   geom_point() + labs(x = "Gene", y = "Cancer Types") + 
#   scale_size_manual(values = sizeLable, guide = "legend", name = "Gastic q value") + 
#   scale_fill_manual(values = fillLabel, guide = "legend", name = "SCNA levels") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



##############

nrStat <- reshape2::melt(nrPeakStat, measure.vars = c('Amplification', 'Deletion')) %>% 
  mutate(Approved.symbol = factor(Approved.symbol, levels = unique(nrPeakQvalue$Approved.symbol)))

plot2 <- ggplot(nrStat, aes(x = Approved.symbol, y = value, fill = variable)) +
  geom_bar(stat = "identity")  + labs(y = 'Number of dysregulated genes') + 
  scale_y_continuous(position = "right") + 
  scale_fill_manual(values = c('Deletion' = 'blue', 'Amplification' = "red"), name = '') + 
  theme(legend.position = "top", legend.justification = "left", 
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.title.y = element_blank()) + coord_flip()


# plots2 <- ggplot(subset(nrStat, numOfDis > 5), aes(x = reorder(Approved.symbol, -index), y = value, fill = variable)) +
#   geom_bar(stat = "identity")  + labs(y = 'Number of dysregulated genes') + 
#   scale_y_continuous(position = "right") + 
#   scale_fill_manual(values = c('Deletion' = 'blue', 'Amplification' = "red"), name = '') + 
#   theme(legend.position = "top", legend.justification = "left", 
#         panel.background = element_rect(fill = "white", color = "black"), 
#         axis.title.y = element_blank()) + coord_flip()


##############
enrichResideInPeak <- readRDS(file = '/data/enrichResideInPeak.rds')
enrichResideInPeak <- enrichResideInPeak %>% mutate(OR = ifelse(Direction == 'Deletion', -OR, OR))


plot3 <- ggplot(enrichResideInPeak, aes(x = Disease, y = OR, fill = Direction)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") + theme_minimal() + 
  scale_fill_manual(values = c('Deletion' = 'blue', 'Amplification' = "red"), name = '') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "top", legend.justification = "left") + 
  geom_text(data = subset(enrichResideInPeak, pValue < 0.05 & abs(OR) > 1),
            size = 3, check_overlap = TRUE, aes(label = pValue), nudge_y = 0.005)


pdf(file = '/result/section2/nrGenePeakEnrich_1.pdf')

plot1
plot2
plot3

dev.off()


# pdf(file = '/result/section2/nrGenePeakEnrichSupple_1.pdf')
# 
# plots1
# plots2
# 
# dev.off()
# 



