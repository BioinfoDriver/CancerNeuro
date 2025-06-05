library(dplyr)
library(DESeq2)
library(tibble)
library('GEOquery')


############################################# LIHC, ICGC
icgc.rnaseq.sample <- readRDS('/data/curatedRNAseq/LIHC/ICGC-HCC/icgc_rnaseq_filter_sample.rds')
load(file='/data/curatedRNAseq/LIHC/ICGC-HCC/icgc_lihc_filter_rnaseq.RData')

#replace all NA values with zero
icgc.lihc.filter.count <- icgc.lihc.filter.count %>% replace(is.na(.), 0)


# sample information
icgc.rnaseq.sample$submitted_sample_id <- do.call(rbind, strsplit(icgc.rnaseq.sample$submitted_sample_id, split='_'))[, 1]
tumor.sams <- subset(icgc.rnaseq.sample, specimen_type == 'Tumor')$submitted_sample_id
normal.sams <- subset(icgc.rnaseq.sample, specimen_type == 'Normal')$submitted_sample_id
paired.sams <- setdiff(intersect(tumor.sams, normal.sams), c('RK062', 'RK080')) # 173
# paired.sams <- intersect(tumor.sams, normal.sams)

# primary and normal samples
tumor.sams <- subset(icgc.rnaseq.sample, submitted_sample_id %in% paired.sams & specimen_type == 'Tumor')
normal.sams <- subset(icgc.rnaseq.sample, submitted_sample_id %in% paired.sams & specimen_type == 'Normal')
paired.sams.count <- icgc.lihc.filter.count[, c(tumor.sams$icgc_sample_id, normal.sams$icgc_sample_id)]


# sample information
paired.sams.info <- rbind(tumor.sams, normal.sams) %>% select(icgc_sample_id, submitted_sample_id, specimen_type)
paired.sams.info <- paired.sams.info %>% remove_rownames() %>% column_to_rownames(var = 'icgc_sample_id')


# > dim(paired.sams.count)
# [1] 22913   346
# Filtering
icgc.exp.genes <- rownames(paired.sams.count)[rowSums(paired.sams.count >= 1) >= 70] 
paired.sams.count <- paired.sams.count[icgc.exp.genes, ]
# > dim(paired.sams.count)
# [1] 20077   346

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = paired.sams.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)


save(paired.sams.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/LIHC_ICGC.RData')



############################################# THCA, GSE213647
sam.info <- readxl::read_xlsx(path = '/data/curatedRNAseq/THCA/GSE213647/41467_2024_45366_MOESM5_ESM.xlsx', sheet = 13, col_names = T)
sam.info$PatientID <- sub("-.*", "", sam.info$`Patient ID`)
sam.info <- sam.info %>% subset(`Library kit` == 'TruSeq RNA access library') %>% dplyr::count(PatientID) %>% subset(n == 2) %>% select(PatientID) %>% left_join(sam.info, by = 'PatientID')

stat.info <- sam.info %>% group_by(PatientID) %>% summarise(ntype = n_distinct(`Tissue type`), nkit = n_distinct(`Library kit`), 
                                                      nage = n_distinct(Age), nyear = n_distinct(Dx_year), nsex = n_distinct(Sex))
# 434
sam.info <- stat.info %>% subset(ntype == 2) %>% select(PatientID) %>% left_join(sam.info, by = 'PatientID') 

# count data
f.names <- list.files(path = "/data/curatedRNAseq/THCA/GSE213647/GSE213647_RAW/", 
                 pattern = "ReadsPerGene.out.tab", full.names = TRUE)

exp.count <- lapply(f.names, function(f) read.table(file = f, skip = 4))

exp <- as.data.frame(sapply(exp.count, function(x) x[, 2]))
rownames(exp) <- exp.count[[1]][, 1]
colnames(exp) <- sub("_.*", "", gsub(pattern = '/data/curatedRNAseq/THCA/GSE213647/GSE213647_RAW/', replacement = '', f.names))

exp <- exp[, sam.info$`GEO ID`]


# batch effect corrrction using ComBat_seq
batch <- sam.info %>% select("GEO ID", "Fixation type") %>% tibble::column_to_rownames(var = "GEO ID")
exp <- as.matrix(exp)
bc.exp.count <- sva::ComBat_seq(counts = exp, batch=batch$`Fixation type`, group=NULL)


# > dim(bc.exp.count)
# [1] 58288   104

# Filtering
exp.genes <- rownames(bc.exp.count)[rowSums(bc.exp.count >= 1) >= 21] 
bc.exp.count <- bc.exp.count[exp.genes, ]

# > dim(bc.exp.count)
# [1] 35015   104

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = bc.exp.count, colData=batch, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
bc.exp.norm <- assay(vsd)

save(bc.exp.count, bc.exp.norm, sam.info, file = '/data/curatedRNAseq/THCA_GSE213647.RData')


############################################# STAD, GSE179252
gse <- getGEO(filename = '/data/curatedRNAseq/STAD/GSE179252/GSE179252_series_matrix.txt.gz')
cli.data <- pData(gse) 

cli.data <- cli.data[, c('title', 'geo_accession', "age:ch1", "gender:ch1", "tissue:ch1", "tumor stage:ch1")]
colnames(cli.data) <- c('title', 'geo_accession', 'age', "gender", "tissue.type", "tumor.stage")

tumor.sams <- subset(cli.data, tissue.type == 'gastric tumor') %>% mutate(PatientID = paste0('P', seq(length(title))))
norm.sams <- subset(cli.data, tissue.type == 'paired normal gastric tissue') %>% mutate(PatientID = paste0('P', seq(length(title))))

paired.sams.info <- rbind(tumor.sams, norm.sams)
rownames(paired.sams.info) <- paired.sams.info$title

exp.count <- read.csv(file = '/data/curatedRNAseq/STAD/GSE179252/GSE179252_Raw_gene_counts_matrix.txt.gz', sep = '\t', header = T)
exp.count <- exp.count %>% column_to_rownames(var = 'gene_id')


# > dim(exp.count)
# [1] 58735    76
# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 15] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 32271    76

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)


save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/STAD_GSE179252.RData')


############################################# ESCA, GSE234304
gse <- getGEO(filename = '/data/curatedRNAseq/ESCA/GSE234304/GSE234304_series_matrix.txt.gz')
cli.data <- pData(gse) 

cli.data <- cli.data[, c('title', 'geo_accession', "tissue:ch1")]
colnames(cli.data) <- c('title', 'geo_accession', "tissue.type")
cli.data <- subset(cli.data, tissue.type != 'Barretts esophagus')

cli.data <- cli.data %>% mutate(PatientID = paste0('P', gsub(pattern = 'N|T', '', title)))

cli.data <- cli.data %>% group_by(PatientID) %>% dplyr::count(PatientID) %>% subset(n == 2) %>% select(PatientID) %>% left_join(cli.data, by = 'PatientID')

# count
exp.count <- read.csv(file = '/data/curatedRNAseq/ESCA/GSE234304/GSE234304_raw_counts.csv', header = TRUE, quote = '\t')
exp.count <- exp.count %>% mutate(X.ID. = gsub('"', '', X.ID.)) %>% column_to_rownames(var = 'X.ID.') %>% select(-X.Gene.name.)

colnames(exp.count) <- gsub("^X\\.X|\\.$", "", colnames(exp.count))
exp.count <- exp.count[, cli.data$title]


# > dim(exp.count)
# [1] 57500    20

# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 4] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 39333    20

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=cli.data, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)


save(exp.count, paired.sams.norm, cli.data, file = '/data/curatedRNAseq/ESCA_GSE234304.RData')


############################################# KIRC, GSE151419
exp.count <- read.csv(file = '/data/curatedRNAseq/KIRC/GSE151419/GSE151419_RNA-Seq_RAW_counts.txt', header = T, sep = '\t')
exp.count <- exp.count %>% column_to_rownames(var = 'Geneid') %>% select(-Chr, -Start, -End, -Strand, -Length)

colnames(exp.count) <- gsub(pattern = '^X', '', colnames(exp.count))

tumor.sams <- grep('C', colnames(exp.count), value = T)
tumor.sams <- gsub('C', '', tumor.sams)

norm.sams <- grep('N', colnames(exp.count), value = T)
norm.sams <- gsub('N', '', norm.sams)

paired.sams <- intersect(norm.sams, tumor.sams)

paired.sams.info <- data.frame(PatientID = c(paired.sams, paired.sams), 
                               sampleID = c(paste0(paired.sams, 'C'), paste0(paired.sams, 'N')), sample.type = rep(c('Tumor', 'Normal'), each = length(paired.sams)))

exp.count <- exp.count[, paired.sams.info$sampleID]

# > dim(exp.count)
# [1] 58735    28

# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 6] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 27601    28

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)


save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/KIRC_GSE151419.RData')


############################################# KIRP, GSE180777
exp.count <- read.csv(file = '/data/curatedRNAseq/KIRP/GSE180777/GSE180777_processed_data_files.csv', header = T, sep = ',')

exp.count <- exp.count %>% mutate(gene_id = sub("\\|.*", "", gene_id)) %>% column_to_rownames(var = 'gene_id')


colnames(exp.count) <- gsub(pattern = '^X', '', colnames(exp.count))

tumor.sams <- grep('T', colnames(exp.count), value = T)
tumor.sams <- gsub('T', '', tumor.sams)

norm.sams <- grep('N', colnames(exp.count), value = T)
norm.sams <- gsub('N', '', norm.sams)

paired.sams <- intersect(norm.sams, tumor.sams)

paired.sams.info <- data.frame(PatientID = c(paired.sams, paired.sams), 
                               sampleID = c(paste0(paired.sams, 'T'), paste0(paired.sams, 'N')), sample.type = rep(c('Tumor', 'Normal'), each = length(paired.sams)))

exp.count <- exp.count[, paired.sams.info$sampleID]


# > dim(exp.count)
# [1] 60675   106

# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 21] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 50172   106

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)


save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/KIRP_GSE180777.RData')



############################################# PRAD, GSE246067
exp.count <- readxl::read_excel(path = "/data/curatedRNAseq/PRAD/GSE246067/GSE246067_ForGEO_RNAseq_Counts_115cases.xlsx", col_names = T, trim_ws = T)

exp.count <- exp.count %>% column_to_rownames(var = 'gene_id') %>% select(-refseq_id)

paired.sams.info <- data.frame(sampleID = colnames(exp.count )) %>% 
  mutate(PatientID = stringr::str_split_i(sampleID, pattern = '_', i = 1), tissue.type = stringr::str_split_i(sampleID, pattern = '_', i = 2))

paired.sams.info <- paired.sams.info %>% column_to_rownames(var = 'sampleID')


exp.count <- exp.count[, rownames(paired.sams.info)]

# > dim(exp.count)
# [1] 26467   230

# Filtering
exp.count <- round(exp.count)
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 46] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 20260   230

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)

save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/PRAD_GSE246067.RData')


############################################# HNSC, GSE178537
exp.count <- read.csv(file = '/data/curatedRNAseq/HNSC/GSE178537/GSE178537_HNSCC_Patient_expected_count.txt', header = T, sep = '\t')

exp.count <- exp.count %>% mutate(gene_id = stringr::str_split_i(pattern = '_', i = 1, gene_id))
exp.count <- exp.count %>% column_to_rownames(var = 'gene_id')
exp.count <- round(exp.count)


paired.sams.info <- data.frame(sampleID = colnames(exp.count)) %>% mutate(tissue.type = stringr::str_split_i(pattern = '\\.', i = 2, sampleID), 
                                                                          PatientID = gsub(pattern = 'X', '', stringr::str_split_i(pattern = '\\.', i = 1, sampleID)))

paired.sams.info <- paired.sams.info %>% mutate(PatientID = paste0('P', PatientID))
paired.sams.info <- subset(paired.sams.info, tissue.type %in% c('IT', 'N'))

paired.sams.info <- paired.sams.info %>% group_by(PatientID) %>% dplyr::count(PatientID) %>% subset(n == 2) %>% select(PatientID) %>% left_join(paired.sams.info, by = 'PatientID')
paired.sams.info <- paired.sams.info %>% column_to_rownames(var = 'sampleID')

exp.count <- exp.count[, rownames(paired.sams.info)]
exp.count <- round(exp.count)

# > dim(exp.count)
# [1] 58676    38

# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 8] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 33799    38

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)


save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/HNSC_GSE178537.RData')


############################################# UCEC, GSE183185
exp.count <- readxl::read_excel(path = '/data/curatedRNAseq/UCEC/GSE183185/GSE183185_counts_anno.xls', col_names = T, trim_ws = T)
exp.count <- exp.count %>% column_to_rownames(var = 'gene_id') %>% select(-Dbxref, -product, -GO_id, -GO_term, -pathway, -pathway_description, -X)


paired.sams.info <- data.frame(sampleID = colnames(exp.count)) %>% mutate(PatientID = paste0('P', stringr::str_replace(sampleID, "^.", "")), 
                                                                          tissue.type = substr(sampleID, 1, 1))

paired.sams.info <- paired.sams.info %>% column_to_rownames(var = 'sampleID') %>% mutate(tissue.type = ifelse(tissue.type == 'C', 'Tumor', 'Normal'))
exp.count <- exp.count[, rownames(paired.sams.info)]

# > dim(exp.count)
# [1] 20030    60

# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 12] 
exp.count <- exp.count[exp.genes, ]
exp.count <- round(exp.count)

# > dim(exp.count)
# [1] 18008    60

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)

save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/UCEC_GSE183185.RData')



############################################# COAD, Nunes_Nature_2024
cli.data <- readxl::read_excel(path = '/data/curatedRNAseq/CRCdata/Nunes_Nature_2024/Supplementary_Table_01.xlsx', skip = 2, col_names = T, trim_ws = T)
cli.data <- cli.data %>% subset(`Primary Site Disease` == 'Colon' & `Histology Subtype` == 'Adenocarcinoma' & `Pre-Treated` == 'Untreated')
cli.data <- cli.data %>% subset(`RNA Normal Sample Barcode` != 'Not_Applicable')

paired.sams.info <- data.frame(PatientID = rep(cli.data$`Sample ID`, 2), sampleID = c(cli.data$`RNA Tumor Sample Barcode`, cli.data$`RNA Normal Sample Barcode`), 
                               tissue.type = rep(c('Tumor', 'Normal'), each = nrow(cli.data)))
paired.sams.info <- paired.sams.info %>% column_to_rownames(var = 'sampleID')

# exp
exp.count <- read.csv(file = '/data/curatedRNAseq/CRCdata/Nunes_Nature_2024/CRC.SW.mRNA.count.txt', header = T, sep = '\t')
exp.count <- exp.count %>% column_to_rownames(var = 'ENSG.ID')
exp.count <- exp.count[, rownames(paired.sams.info)]

# > dim(exp.count)
# [1] 19765   122

# Filtering
exp.genes <- rownames(exp.count)[rowSums(exp.count >= 1) >= 24] 
exp.count <- exp.count[exp.genes, ]

# > dim(exp.count)
# [1] 17835   122

# The variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=paired.sams.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
paired.sams.norm <- assay(vsd)

save(exp.count, paired.sams.norm, paired.sams.info, file = '/data/curatedRNAseq/COAD_Nunes_Nature_2024.RData')



############################################# BLCA, GSE229410
# Normalized (RPKM)
norm.exp <- read.csv(file = '/data/curatedRNAseq/BLCA/GSE229410/GSE229410_Processed_data.txt', header = T, sep = '\t', fileEncoding = 'UTF-16')
norm.exp <- norm.exp %>% column_to_rownames(var = 'SYMBOL1')


paired.sams.info <- data.frame(sampleID = colnames(norm.exp)) %>% mutate(PatientID = gsub(pattern = 'A|P', replacement = '', sampleID))
paired.sams.info <- paired.sams.info %>% group_by(PatientID) %>% dplyr::count(PatientID) %>% subset(n == 2) %>% select(PatientID) %>% left_join(paired.sams.info, by = 'PatientID')

paired.sams.info$tissue.type <- 'Tumor'
paired.sams.info$tissue.type[grep(pattern = 'P', paired.sams.info$sampleID)] <- 'Normal'


norm.exp <- norm.exp[, paired.sams.info$sampleID]
paired.sams.info <- paired.sams.info %>% column_to_rownames('sampleID')

# > dim(norm.exp)
# [1] 19664    68

# Filtering
norm.exp <- norm.exp[rowSums(norm.exp >= 0.1) >= 14, ]
norm.exp <- log2(norm.exp + 1)

# > dim(norm.exp)
# [1] 15407    68

save(norm.exp, paired.sams.info, file = '/data/curatedRNAseq/BLCA_GSE229410.RData')


############################################# BRCA, GSE233242
# Normalized (TPM)
gse <- getGEO(filename = '/data/curatedRNAseq/BRCA/GSE233242/GSE233242_series_matrix.txt')
cli.data <- pData(gse) 

cli.data <- cli.data[, c('title', 'geo_accession', "tissue:ch1", "Sex:ch1", "subtype:ch1", "description")]
colnames(cli.data) <- c('title', 'geo_accession', "tissue.type", "sex", 'subtype', "sampleID")

cli.data <- cli.data %>% mutate(PatientID = paste0('P', stringr::str_split_i(pattern = ' ', title, i = 2)))
cli.data <- cli.data %>% remove_rownames() %>% column_to_rownames(var = 'sampleID')


# exp
norm.exp <- read.csv(file = "/data/curatedRNAseq/BRCA/GSE233242/GSE233242_Processed_file_tpm.csv", header = T, sep = ',')
norm.exp <- norm.exp %>% column_to_rownames('gene_id') %>% select(-gene_name, -transcriptLen)

# > dim(norm.exp)
# [1] 15044    86

# Filtering
norm.exp <- norm.exp[rowSums(norm.exp >= 0.1) >= 17, ]
norm.exp <- log2(norm.exp + 1)

# > dim(norm.exp)
# [1] 15042    86

colnames(norm.exp) <- rownames(cli.data)
save(norm.exp, cli.data, file = '/data/curatedRNAseq/BRCA_GSE233242.RData')



############################################# KIRC, GSE251905
# Normalized (log2(CPM))
clinical.data <- readxl::read_excel(path = "/data/curatedRNAseq/KIRC/GSE251905/GSE251905_Extended_Data_Table_2.xlsx", 
                                    skip = 1, col_names = F, n_max = 3, trim_ws = T)
colnames(clinical.data) <- paste0('V', seq(ncol(clinical.data)))

clinical.data <- clinical.data %>% column_to_rownames(var = 'V1')
clinical.data <- t(clinical.data) %>% as.data.frame()


clinical.data <- subset(clinical.data, `Tissue Type` != 'Metastasis' & !(Diagnosis %in% c('Chromophobe', 'FH Def RCC', 'Oncocytic', 'Papillary')))
colnames(clinical.data) <- c('PatientID', 'TissueType', 'Diagnosis')
clinical.data <- clinical.data %>% rownames_to_column(var = 'sampleID')

paired.sams.info <- clinical.data %>% group_by(PatientID) %>% dplyr::count(PatientID) %>% subset(n == 2) %>% select(PatientID) %>% left_join(clinical.data, by = 'PatientID') 
paired.sams.info <- paired.sams.info %>% column_to_rownames(var = 'sampleID')

# exp data
norm.exp <- readxl::read_excel(path = "/data/curatedRNAseq/KIRC/GSE251905/GSE251905_Extended_Data_Table_2.xlsx", 
                               skip = 4, col_names = F, trim_ws = T)
colnames(norm.exp) <- paste0('V', seq(ncol(norm.exp)))
norm.exp <- norm.exp %>% column_to_rownames(var = 'V1')

norm.exp <- norm.exp[, rownames(paired.sams.info)]

# > dim(norm.exp)
# [1] 24438    54

# Filtering
norm.exp <- norm.exp[rowSums(norm.exp >= 1) >= 11, ]

# > dim(norm.exp)
# [1] 22352    54

save(norm.exp, paired.sams.info, file = '/data/curatedRNAseq/KIRC_GSE251905.RData')


############################################# LUSC, Shankha_CELL_2021
# Log2-transformed upper-quartile (UQ)-normalized FPKM values
# Genes (RNA-seq), present in fewer than 30% of samples (i.e., missing in > 70% of samples) were removed from the respective datasets.
norm.exp <- readxl::read_excel(path = '/data/curatedRNAseq/LUSC/Shankha_CELL_2021/mmc2.xlsx', col_names = T, trim_ws = T, sheet = 5, na = 'NA')
norm.exp <- norm.exp %>% column_to_rownames('ENSEMBL') %>% select(-id, -gene_id)

# 108 tumor samples were used in the analysis
cli.data <- readxl::read_excel(path = '/data/curatedRNAseq/LUSC/Shankha_CELL_2021/mmc1.xlsx', col_names = T, trim_ws = T, sheet = 2)
cli.data <- cli.data %>% select(Sample.ID, Participant, Type)

cli.data <- subset(cli.data, Sample.ID %in% colnames(norm.exp))
cli.data <- cli.data %>% group_by(Participant) %>% dplyr::count(Participant) %>% subset(n == 2) %>% select(Participant) %>% left_join(cli.data, by = 'Participant')
colnames(cli.data) <- c('PatientID', 'Sample.ID', 'tissue.type')
cli.data <- cli.data %>% column_to_rownames('Sample.ID')

norm.exp <- norm.exp[, rownames(cli.data)]
norm.exp[is.na(norm.exp)] <- 0

# > dim(norm.exp)
# [1] 21792   188

# Filtering
norm.exp <- norm.exp[rowSums(norm.exp >= 1) >= 38, ]
# > dim(norm.exp)
# [1] 21792   188

save(norm.exp, cli.data, file = '/data/curatedRNAseq/LUSC_Shankha_CELL_2021.RData')


############################################# LUAD, GSE40419
# Normalized (RPKM)
exp.norm <- read.csv(file = '/data/curatedRNAseq/LUAD/GSE40419/GSE40419_LC-87_RPKM_expression.txt', sep = '\t', header = T)

norm.sams <- grep(pattern = '_nor', colnames(exp.norm)[-c(1:6)], value = T)
tumor.sams <- setdiff(colnames(exp.norm)[-c(1:6)], norm.sams)
norm.sams <- gsub(pattern = '_nor', '', norm.sams)

norm.sams <- intersect(norm.sams, tumor.sams)
tumor.sams <- intersect(norm.sams, tumor.sams)


# sample information
paired.sams.info <- data.frame(labels=c(rep('Tumor', length(tumor.sams)), rep('Normal', length(norm.sams))), PatientID = c(tumor.sams, norm.sams))
rownames(paired.sams.info) <- c(tumor.sams, paste0(norm.sams, '_nor'))


exp.norm <- exp.norm[, c(colnames(exp.norm)[1:6], rownames(paired.sams.info))]

# When multiple probes corresponded to the same gene, 
# the probe with the highest normalized intensity averaged over all samples was used.
exp.norm <- lapply(split.data.frame(exp.norm[, -c(2:6)], f = ~ gene), function(x){
  
  x <- x[which.max(rowSums(x[, -1])), ]
  return(x)
})

exp.norm <- do.call(rbind, exp.norm) 
exp.norm <- exp.norm %>% select(-gene)

# > dim(exp.norm)
# [1] 22427   154

# Filtering
exp.norm <- exp.norm[rowSums(exp.norm >= 0.1) >= 31, ]
exp.norm <- log2(exp.norm + 1)
# > dim(exp.norm)
# [1] 17863   154

save(exp.norm, paired.sams.info, file = '/data/curatedRNAseq/LUAD_GSE40419.RData')

