
library('dplyr')
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

#################
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")

anno450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- anno450k %>% as.data.frame() %>% subset(UCSC_RefGene_Name != '')


splitData <- function(data){
  
  RefGene_Name <- unlist(strsplit(data$UCSC_RefGene_Name, split = ';'))
  RefGene_Accession <- unlist(strsplit(data$UCSC_RefGene_Accession, split = ';'))
  RefGene_Group  <- unlist(strsplit(data$UCSC_RefGene_Group, split = ';'))
  
  
  data <- data.frame(Name = data$Name[1], RefGene_Name, RefGene_Accession, RefGene_Group)
  
  return(data) 
}

anno450kGeneAnno <- anno450k %>% group_by(Name) %>% do(splitData(.))

anno450k <- anno450k %>% left_join(anno450kGeneAnno, by = 'Name')

saveRDS(anno450k, file = '/data/methAnno450k.rds')
