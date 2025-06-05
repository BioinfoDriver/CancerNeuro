library(dplyr)
######################
selectRes <- readRDS(file = '/data/nRSelectScoreRes.rds')

durgTargetGenes <- read.csv(file = '/data/oncokb_biomarker_drug_associations.tsv', header = T, sep = '\t', stringsAsFactors = F)
geneInfo <- read.csv(file= '/data/gene_with_protein_product.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


durgTargetGenes <- durgTargetGenes %>% distinct(Gene, .keep_all = TRUE)

sigselectRes <- selectRes %>% subset(FDR & APC_good) %>% mutate(SFE_1 = gsub("_.*$", "", SFE_1), SFE_2 = gsub("_.*$", "", SFE_2)) %>% 
  subset(SFE_1 %in% durgTargetGenes$Gene | SFE_2 %in% durgTargetGenes$Gene) %>% subset(!(SFE_1 %in% durgTargetGenes$Gene & SFE_2 %in% durgTargetGenes$Gene))


sigselectRes <- sigselectRes %>% inner_join(geneInfo[, c('symbol', 'location')], by = join_by(SFE_1 == symbol)) %>% 
  rename(SFE_1_location = location) %>% inner_join(geneInfo[, c('symbol', 'location')], by = join_by(SFE_2 == symbol)) %>% 
  rename(SFE_2_location = location) 


# > dim(sigselectRes)
# [1] 442  28

######################
# Filtering
sigselectRes <- sigselectRes %>% subset(!((SFE_1_location == SFE_2_location) & 
                                            (!int_type %in% c('AMP - DEL', 'AMP - MUT', 'MUT - MUT'))))
# > dim(sigselectRes)
# [1] 423  30


# Filtering
sigselectRes <- sigselectRes %>% mutate(SFE_1_location_1 = gsub("\\..*$", "", SFE_1_location), SFE_2_location_1 = gsub("\\..*$", "", SFE_2_location))

sigselectRes <- sigselectRes %>% subset(!((SFE_1_location_1 == SFE_2_location_1) & 
                                            (!int_type %in% c('AMP - DEL', 'AMP - MUT', 'MUT - MUT'))))
# > dim(sigselectRes)
# [1] 376  30


# Filtering
sigselectRes <- sigselectRes %>% mutate(SFE_1_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_1_location_1),
                                        SFE_2_location_2 = gsub("([a-zA-Z]+)([0-9]+)", "\\1", SFE_2_location_1),
                                        SFE_1_location_2_subloca = as.numeric(stringr::str_split(SFE_1_location_1, pattern = 'p|q', simplify = T)[, 2]), 
                                        SFE_2_location_2_subloca = as.numeric(stringr::str_split(SFE_2_location_1, pattern = 'p|q', simplify = T)[, 2]))

sigselectRes <- sigselectRes %>% subset(!((SFE_1_location_2 == SFE_2_location_2) & 
                                            (abs(SFE_2_location_2_subloca - SFE_1_location_2_subloca) < 6) & 
                                            (!int_type %in% c('AMP - DEL', 'AMP - MUT', 'MUT - MUT'))))
# > dim(sigselectRes)
# [1] 219  34


# write.table(x = sigselectRes, row.names = F, col.names = T, sep = '\t', quote = F, file = '/result/section2/nrDrugTagetMe.txt')

######################
library(ComplexHeatmap)


showSigselectRes <- sigselectRes %>% subset(freq_1 > 0.040 | freq_2 > 0.012)
meGenePairs <- stringr::str_split(showSigselectRes$name, pattern = ' - ', simplify = TRUE) %>% as.data.frame() %>% rename(Gene1 = V1, Gene2 = V2)
showSigselectRes <- cbind.data.frame(meGenePairs, showSigselectRes)

showSigselectRes <- showSigselectRes %>% mutate(Gene_1 = ifelse(SFE_1 %in% durgTargetGenes$Gene, Gene1, Gene2), 
                                                Gene_2 = ifelse(SFE_1 %in% durgTargetGenes$Gene, Gene2, Gene1)) %>% 
  mutate(APC = ifelse(direction == 'CO', -APC, APC)) %>% distinct(SFE_2_location, .keep_all = T) 
# %>% distinct(SFE_1_location, .keep_all = T) 



showData <- showSigselectRes %>% tidyr::pivot_wider(id_cols = Gene_1, names_from = Gene_2, values_from = APC) %>% 
  column_to_rownames(var = 'Gene_1')



colors<- circlize::colorRamp2(breaks = c(-0.3, 0, 0.006), colors = c("#8BB77B", "white","#80539A"))

p1 <- Heatmap(t(showData), cluster_rows = FALSE, cluster_columns = FALSE, na_col = "white", 
              show_heatmap_legend = F, col = colors,)

lgd <- list(Legend(title = "SELECT score", 
                   col_fun = colors, 
                   at = c(-0.2, -0.15, 0, 0.003, 0.006), 
                   direction = "horizontal"))

pdf(file = '/result/section2/nrDrugTagetMe.pdf')
draw(p1, annotation_legend_list = lgd, annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()


