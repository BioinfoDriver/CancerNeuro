library(limma)
library(tibble)
library(dplyr)


nrs <- readRDS(file = '/data/neurotransmitterReceptors.rds')

# Function
PairedSampleDiffExp <- function(disease, expData, samInfo){
  
  colnames(samInfo) <- c('patientID', 'sampleType')
  
  samInfo$patientID <- factor(samInfo$patientID)
  samInfo$sampleType <- factor(samInfo$sampleType, levels = c('Normal', 'Tumor'))
  
  
  design <- model.matrix(~patientID+sampleType, data = samInfo)
  fit <- lmFit(expData, design)
  fit <- eBayes(fit)
  diffExp <- topTable(fit, coef="sampleTypeTumor", adjust.method = "BH", n = Inf)
  
  diffExp <- diffExp %>% mutate(DISEASE = disease) %>% rownames_to_column(var = "geneSymbol")
  return(diffExp)  
}


setwd('/data/curatedRNAseq')
# grep('.RData', dir(path = '/data/curatedRNAseq'), value = TRUE)
# [1] "BLCA_GSE229410.RData"         "BRCA_GSE233242.RData"         "COAD_Nunes_Nature_2024.RData" "ESCA_GSE234304.RData"         "HNSC_GSE178537.RData"         "KIRC_GSE151419.RData"        
# [7] "KIRC_GSE251905.RData"         "KIRP_GSE180777.RData"         "LIHC_ICGC.RData"              "LUAD_GSE40419.RData"          "LUSC_Shankha_CELL_2021.RData" "PRAD_GSE246067.RData"        
# [13] "STAD_GSE179252.RData"         "THCA_GSE213647.RData"         "UCEC_GSE183185.RData" 
# 

# COAD_Nunes_Nature_2024
load(file = 'COAD_Nunes_Nature_2024.RData')

paired.sams.norm <- paired.sams.norm[intersect(nrs$Ensembl.gene.ID, rownames(paired.sams.norm)), ]
rownames(paired.sams.norm) <- nrs$Approved.symbol[match(rownames(paired.sams.norm), nrs$Ensembl.gene.ID)]

coadDiffExp <- PairedSampleDiffExp('COAD', paired.sams.norm, paired.sams.info)


# ESCA_GSE234304.RData
load(file = 'ESCA_GSE234304.RData')
cli.data <- cli.data %>% column_to_rownames(var = 'title') %>% select(PatientID, tissue.type) %>% mutate(tissue.type = ifelse(tissue.type == 'Tumor', 'Tumor', 'Normal'))

paired.sams.norm <- paired.sams.norm[intersect(nrs$Ensembl.gene.ID, rownames(paired.sams.norm)), ]
rownames(paired.sams.norm) <- nrs$Approved.symbol[match(rownames(paired.sams.norm), nrs$Ensembl.gene.ID)]

escaDiffExp <- PairedSampleDiffExp('ESCA', paired.sams.norm, cli.data)


# HNSC_GSE178537.RData
load(file = 'HNSC_GSE178537.RData')

paired.sams.info <- paired.sams.info %>% mutate(tissue.type = ifelse(tissue.type == 'N', 'Normal', 'Tumor'))

paired.sams.norm <- paired.sams.norm[intersect(nrs$Ensembl.gene.ID, rownames(paired.sams.norm)), ]
rownames(paired.sams.norm) <- nrs$Approved.symbol[match(rownames(paired.sams.norm), nrs$Ensembl.gene.ID)]

hnscDiffExp <- PairedSampleDiffExp('HNSC', paired.sams.norm, paired.sams.info)


# KIRP_GSE180777.RData
load(file = 'KIRP_GSE180777.RData')

paired.sams.info <- paired.sams.info %>% column_to_rownames(var = 'sampleID')

paired.sams.norm <- paired.sams.norm[intersect(nrs$Ensembl.gene.ID, rownames(paired.sams.norm)), ]
rownames(paired.sams.norm) <- nrs$Approved.symbol[match(rownames(paired.sams.norm), nrs$Ensembl.gene.ID)]

kirpDiffExp <- PairedSampleDiffExp('KIRP', paired.sams.norm, paired.sams.info)


# LIHC_ICGC.RData
load(file = 'LIHC_ICGC.RData')

paired.sams.norm <- paired.sams.norm[intersect(nrs$Approved.symbol, rownames(paired.sams.norm)), ]

lihcDiffExp <- PairedSampleDiffExp('LIHC', paired.sams.norm, paired.sams.info)


# PRAD_GSE246067.RData
load(file = 'PRAD_GSE246067.RData')
paired.sams.info <- paired.sams.info %>% mutate(tissue.type = ifelse(tissue.type == 'Benign', 'Normal', 'Tumor'))

paired.sams.norm <- paired.sams.norm[intersect(nrs$Approved.symbol, rownames(paired.sams.norm)), ]

pradDiffExp <- PairedSampleDiffExp('PRAD', paired.sams.norm, paired.sams.info)


# STAD_GSE179252.RData
load(file = 'STAD_GSE179252.RData')
paired.sams.info <- paired.sams.info %>% select(PatientID, tissue.type) %>% mutate(tissue.type = ifelse(tissue.type == 'gastric tumor', 'Tumor', 'Normal'))

paired.sams.norm <- paired.sams.norm[intersect(nrs$Ensembl.gene.ID, rownames(paired.sams.norm)), ]
rownames(paired.sams.norm) <- nrs$Approved.symbol[match(rownames(paired.sams.norm), nrs$Ensembl.gene.ID)]

stadDiffExp <- PairedSampleDiffExp('STAD', paired.sams.norm, paired.sams.info)



# "THCA_GSE213647.RData"  
load(file = "THCA_GSE213647.RData")

sam.info <- sam.info[, c('PatientID', "GEO ID", "Tissue type")]
colnames(sam.info) <- c('PatientID', "GEO.ID", "Tissue.type")
sam.info <- sam.info %>% column_to_rownames(var = 'GEO.ID') %>% mutate(Tissue.type = ifelse(Tissue.type == 'Normal', 'Normal', 'Tumor'))


rownames(bc.exp.norm) <- stringr::str_split_i(rownames(bc.exp.norm), '\\.', i = 1)
bc.exp.norm <- bc.exp.norm[intersect(nrs$Ensembl.gene.ID, rownames(bc.exp.norm)), ]
rownames(bc.exp.norm) <- nrs$Approved.symbol[match(rownames(bc.exp.norm), nrs$Ensembl.gene.ID)]

thcaDiffExp <- PairedSampleDiffExp('THCA', paired.sams.norm, paired.sams.info)


# "UCEC_GSE183185.RData" 
load(file = "UCEC_GSE183185.RData")

paired.sams.norm <- paired.sams.norm[intersect(nrs$Approved.symbol, rownames(paired.sams.norm)), ]

ucecDiffExp <- PairedSampleDiffExp('UCEC', paired.sams.norm, paired.sams.info)



# "BLCA_GSE229410.RData"
load(file = "BLCA_GSE229410.RData")

norm.exp <- norm.exp[intersect(nrs$Approved.symbol, rownames(norm.exp)), ]
blcaDiffExp <- PairedSampleDiffExp('BLCA', norm.exp, paired.sams.info)


# "BRCA_GSE233242.RData"
load(file = "BRCA_GSE233242.RData")
cli.data <- cli.data %>% select(PatientID, subtype) %>% mutate(subtype = ifelse(subtype == 'Normal', 'Normal', 'Tumor'))

norm.exp <- norm.exp[intersect(nrs$Ensembl.gene.ID, rownames(norm.exp)), ]
rownames(norm.exp) <- nrs$Approved.symbol[match(rownames(norm.exp), nrs$Ensembl.gene.ID)]

brcaDiffExp <- PairedSampleDiffExp('BRCA', norm.exp, cli.data)


# "KIRC_GSE251905.RData"
load(file = "KIRC_GSE251905.RData")
paired.sams.info <- paired.sams.info %>% select(PatientID, TissueType) %>% mutate(TissueType = ifelse(TissueType == 'Primary Tumour', 'Tumor', 'Normal'))

norm.exp <- norm.exp[intersect(nrs$Approved.symbol, rownames(norm.exp)), ]

kircDiffExp <- PairedSampleDiffExp('KIRC', norm.exp, paired.sams.info)


# LUSC_Shankha_CELL_2021.RData
load(file = "LUSC_Shankha_CELL_2021.RData")

cli.data <- cli.data %>% mutate(tissue.type = ifelse(tissue.type == 'Tumor', 'Tumor', 'Normal'))
norm.exp <- norm.exp[intersect(nrs$Ensembl.gene.ID, rownames(norm.exp)), ]
rownames(norm.exp) <- nrs$Approved.symbol[match(rownames(norm.exp), nrs$Ensembl.gene.ID)]

luscDiffExp <- PairedSampleDiffExp('LUSC', norm.exp, cli.data)


# LUAD_GSE40419.RData
load(file = "LUAD_GSE40419.RData")
paired.sams.info <- paired.sams.info %>% select(PatientID, labels)
exp.norm <- exp.norm[intersect(nrs$Approved.symbol, rownames(exp.norm)), ]

luadDiffExp <- PairedSampleDiffExp('LUAD', exp.norm, paired.sams.info)


################################
nRDiffExp <- rbind.data.frame(blcaDiffExp, brcaDiffExp, coadDiffExp, escaDiffExp,           
hnscDiffExp, kircDiffExp, kirpDiffExp, lihcDiffExp, luadDiffExp, luscDiffExp, 
pradDiffExp, stadDiffExp, thcaDiffExp, ucecDiffExp)


saveRDS(nRDiffExp, file = '/data/curatedPanCanNrDiffExp.rds')

