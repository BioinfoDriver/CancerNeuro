

##############
panCanSurPvalue <- readRDS(file = '/data/panCanNrSurPvalue.rds')

# panCanSurPvalue <- subset(panCanSurPvalue, !is.na(HR))
# 
# 
# write.table(panCanSurPvalue, file = '/result/section4/nrReceptorAssociatedWithOS.txt',
#             sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)


############
selectSurPvalue <- subset(panCanSurPvalue, coxPvalue < 0.05 & logrankPvalue < 0.05)
selectSurPvalue <- selectSurPvalue %>% mutate(hrDirect = ifelse(HR > 1.0, 'High level-> worse survival', 'High level-> better survival'))

write.table(selectSurPvalue, file = '/result/section4/nrSigAssociatedWithOS.txt',
             sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)



disDysGCount <- selectSurPvalue %>% group_by(disease, hrDirect) %>% count(name = 'num')
dysGcount <- selectSurPvalue %>% group_by(geneSymbol, hrDirect) %>% count(name = 'num')

dysGcount <- dysGcount %>% group_by(geneSymbol) %>% summarise(totalDysCount = sum(num)) %>% right_join(dysGcount, by = 'geneSymbol')


write.table(disDysGCount, file = '/result/section4/disDysGCount.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)


write.table(dysGcount, file = '/result/section4/dysGcount.txt',
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

