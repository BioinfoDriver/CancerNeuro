
####################TCGA
# CDKN2A/B Deletion, EGFR amplification
gene.cnv.alt <- read.table(file = '/data/all_thresholded.by_genes.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(gene.cnv.alt) <- gene.cnv.alt$Gene.Symbol
gene.cnv.alt <- gene.cnv.alt[, -c(1:3)]
colnames(gene.cnv.alt) <- gsub('\\.', '-', substr(colnames(gene.cnv.alt), 1, 15))
gene.cnv.alt <- as.data.frame(t(gene.cnv.alt))
gene.cnv.alt$CDKN2AB <- ifelse((gene.cnv.alt$CDKN2A == -2) | (gene.cnv.alt$CDKN2B == -2), 1, 0)
gene.cnv.alt$EGFR <- ifelse(gene.cnv.alt$EGFR == 2, 1, 0)


sub.cli.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
tcga_glioma_cli_mol <- readRDS(file = 'D:/CancerNeuroscience/data/GliomaData/tcga_glioma_cli_mol.rds')
tcga_glioma_cli_mol <- tcga_glioma_cli_mol %>% mutate(bcr_patient_barcode = paste0(bcr_patient_barcode, '-01'))
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

sub.cli.data <- merge(sub.cli.data, tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'SUBTYPE')], by.x = 'patient_id', by.y = 'SAMPLE_BARCODE')
sub.cli.data <- merge(sub.cli.data, tcga_glioma_cli_mol[, c('bcr_patient_barcode', 'MGMT_PROMOTER_STATUS')], 
                      by.x = 'patient_id', by.y = 'bcr_patient_barcode')

sub.cli.data$CDKN2AB <- gene.cnv.alt[sub.cli.data$patient_id, 'CDKN2AB']
sub.cli.data$EGFR <- gene.cnv.alt[sub.cli.data$patient_id, 'EGFR']


# > mean(sub.cli.data$age)
# [1] 42.99012
sub.cli.data <- sub.cli.data %>% mutate(age_categ = ifelse(age >= 43, '>= 43', '<43'), gender = factor(gender,levels = c('FEMALE', 'MALE')), 
                                        # race = factor(race, levels = c('BLACK', 'OTHER', 'WHITE')),
                                        histological_type = factor(histological_type, levels = c('Oligodendroglioma', 'Oligoastrocytoma', 'Astrocytoma')), 
                                        histological_grade = factor(histological_grade, levels = c('G2', 'G3')),
                                        mgmt.status = factor(MGMT_PROMOTER_STATUS, levels = c('Unmethylated', 'Methylated')), 
                                        CDKN2AB = factor(CDKN2AB, levels = c(0, 1)), 
                                        EGFR = factor(EGFR, levels = c(0, 1)), 
                                        SUBTYPE = factor(SUBTYPE, levels = c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt')), 
                                        risk.categ = factor(risk.categ , levels = c('Others', 'C4-like')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, age_categ, gender, histological_type, histological_grade, mgmt.status, CDKN2AB, EGFR, SUBTYPE, risk.categ) # race,


# Fisher's Exact Test
p.values <- sapply(colnames(cli.sig.char)[2:9], function(char){
  # print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  # print(table(cli.sig.char[, c(char, 'risk.categ')]))
  # 
  # or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  p.value <- chisq.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
  # return(c(or, p.value))
  return(p.value)
})

names(p.values) <- colnames(cli.sig.char)[2:9]


# age_categ             gender  histological_type histological_grade        mgmt.status            CDKN2AB               EGFR            SUBTYPE 
# 2.096981e-07       6.559272e-01       2.875508e-12       2.907869e-08       1.951244e-21       2.200593e-18       9.669644e-24       8.070818e-56

####Plot
library('ggpubr')

age.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('age_categ', 'risk.categ')])), 
                     "risk.categ", "Freq", fill = "age_categ", color="age_categ", 
                     ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                     lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                            legend.title = element_text(size = 5), legend.text = element_text(size = 5))

grade.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('histological_grade', 'risk.categ')])), 
                      "risk.categ", "Freq", fill = "histological_grade", color="histological_grade", 
                      ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                      lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                             legend.title = element_text(size = 5), legend.text = element_text(size = 5))

histype.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('histological_type', 'risk.categ')])), 
                      "risk.categ", "Freq", fill = "histological_type", color="histological_type", 
                      ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                      lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                             legend.title = element_text(size = 5), legend.text = element_text(size = 5))


mgmt.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('mgmt.status', 'risk.categ')])), 
                         "risk.categ", "Freq", fill = "mgmt.status", color="mgmt.status", 
                         ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                         lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                                legend.title = element_text(size = 5), legend.text = element_text(size = 5))

CDKN2AB.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('CDKN2AB', 'risk.categ')])), 
                   "risk.categ", "Freq", fill = "CDKN2AB", color="CDKN2AB", 
                   ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                   lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                          legend.title = element_text(size = 5), legend.text = element_text(size = 5))

EGFR.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('EGFR', 'risk.categ')])), 
                       "risk.categ", "Freq", fill = "EGFR", color="EGFR", 
                       ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                       lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                              legend.title = element_text(size = 5), legend.text = element_text(size = 5))

SUBTYPE.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('SUBTYPE', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "SUBTYPE", color="SUBTYPE", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))


ggsave(ggarrange(age.p, grade.p, histype.p, mgmt.p, CDKN2AB.p, EGFR.p, SUBTYPE.p, ncol=4, nrow=2), width = 20, height = 12, units = 'cm', 
       file='/result/section5/lgglike/tcga_cli_risk_com.pdf')


####################CGGA 693
sub.cli.data <- readRDS(file='/data/LggRiskScores/plgg_693_risk_score.rds')

# > mean(sub.cli.data$Age, na.rm = T)
# [1] 39.97509
sub.cli.data$Gender <- stringr::str_trim(sub.cli.data$Gender, side = "right")

sub.cli.data <- sub.cli.data %>% mutate(os = Censor..alive.0..dead.1., os_time = OS, 
                                        age_categ = ifelse(Age >= 40, '>= 40', '<40'), gender = factor(Gender,levels = c('Female', 'Male')), 
                                        Histology = factor(Histology, levels = c('O', 'OA', 'AO', 'AOA', 'A', 'AA')), 
                                        Grade = factor(Grade, levels = c('WHO II', 'WHO III')),
                                        idh.status = factor(IDH_mutation_status, levels = c('Mutant', 'Wildtype')), 
                                        X1p19q.codeletion = factor(X1p19q_codeletion_status, levels = c('Codel', 'Non-codel')), 
                                        mgmt.status = factor(MGMTp_methylation_status, levels = c('un-methylated', 'methylated')), 
                                        TMZ.treated = factor(Chemo_status..TMZ.treated.1.un.treated.0., levels = c(0, 1)), 
                                        Radio.status = factor(Radio_status..treated.1.un.treated.0., levels = c(0, 1)), 
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, age_categ, gender, Histology, Grade, idh.status, X1p19q.codeletion, 
                                        mgmt.status, risk.categ)


# Fisher's Exact Test
p.values <- sapply(colnames(cli.sig.char)[2:8], function(char){
  # print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  print(table(cli.sig.char[, c(char, 'risk.categ')]))
  # 
  # or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  p.value <- chisq.test(table(cli.sig.char[, c(char, 'risk.categ')]))$p.value
  # return(c(or, p.value))
  return(p.value)
})

names(p.values) <- colnames(cli.sig.char)[2:8]


# age_categ            gender         Histology             Grade        idh.status X1p19q.codeletion       mgmt.status 
# 7.572598e-01      6.345040e-01      1.201596e-03      1.500184e-01      8.786444e-09      2.874916e-03      5.524186e-01

####Plot
library('ggpubr')

X1p19q.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('X1p19q.codeletion', 'risk.categ')])), 
                   "risk.categ", "Freq", fill = "X1p19q.codeletion", color="X1p19q.codeletion", 
                   ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                   lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                          legend.title = element_text(size = 5), legend.text = element_text(size = 5))

IDH.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('idh.status', 'risk.categ')])), 
                      "risk.categ", "Freq", fill = "idh.status", color="idh.status", 
                      ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                      lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                             legend.title = element_text(size = 5), legend.text = element_text(size = 5))

histype.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('Histology', 'risk.categ')])), 
                        "risk.categ", "Freq", fill = "Histology", color="Histology", 
                        ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                        lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                               legend.title = element_text(size = 5), legend.text = element_text(size = 5))

ggsave(ggarrange(X1p19q.p, IDH.p, histype.p, ncol=4, nrow=1), width = 20, height = 6, units = 'cm', 
       file='/result/section5/lgglike/plgg_693_cli_risk_com.pdf')



####################CGGA 325
sub.cli.data <- readRDS(file='/data/LggRiskScores/plgg_325_risk_score.rds')

# > mean(sub.cli.data$Age, na.rm = T)
# [1] 40.66667

sub.cli.data <- sub.cli.data %>% mutate(os = Censor..alive.0..dead.1., os_time = OS, 
                                        age_categ = ifelse(Age >= 41, '>= 41', '<41'), gender = factor(Gender,levels = c('Female', 'Male')), 
                                        Histology = factor(Histology, levels = c('O', 'OA', 'AO', 'AOA', 'A', 'AA')), 
                                        Grade = factor(Grade, levels = c('WHO II', 'WHO III')),
                                        idh.status = factor(IDH_mutation_status, levels = c('Mutant', 'Wildtype')), 
                                        X1p19q.codeletion = factor(X1p19q_codeletion_status, levels = c('Codel', 'Non-codel')),
                                        mgmt.status = factor(MGMTp_methylation_status, levels = c('un-methylated', 'methylated')), 
                                        TMZ.treated = factor(Chemo_status..TMZ.treated.1.un.treated.0., levels = c(0, 1)),
                                        Radio.status = factor(Radio_status..treated.1.un.treated.0., levels = c(0, 1)),
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, age_categ, gender, Grade, idh.status, X1p19q.codeletion, 
                                        mgmt.status, Histology, risk.categ) # , TMZ.treated, Radio.status


# Fisher's Exact Test
p.values <- sapply(colnames(cli.sig.char)[2:8], function(char){
  # print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  print(table(cli.sig.char[, c(char, 'risk.categ')]))
  # 
  # or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  stat <- table(cli.sig.char[, c(char, 'risk.categ')])
  p.value <- chisq.test(stat[rowSums(stat) > 0, ])$p.value
  # return(c(or, p.value))
  return(p.value)
})

names(p.values) <- colnames(cli.sig.char)[2:8]

# age_categ            gender             Grade        idh.status X1p19q.codeletion       mgmt.status         Histology 
# 2.047859e-02      9.211185e-01      1.229676e-06      1.219451e-09      1.594437e-07      1.913149e-01      2.977885e-10 

####Plot
library('ggpubr')
age.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('age_categ', 'risk.categ')])), 
                      "risk.categ", "Freq", fill = "age_categ", color="age_categ", 
                      ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                      lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

grade.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('Grade', 'risk.categ')])), 
                   "risk.categ", "Freq", fill = "Grade", color="Grade", 
                   ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                   lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                          legend.title = element_text(size = 5), legend.text = element_text(size = 5))

X1p19q.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('X1p19q.codeletion', 'risk.categ')])), 
                      "risk.categ", "Freq", fill = "X1p19q.codeletion", color="X1p19q.codeletion", 
                      ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                      lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                             legend.title = element_text(size = 5), legend.text = element_text(size = 5))

IDH.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('idh.status', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "idh.status", color="idh.status", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))

histype.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('Histology', 'risk.categ')])), 
                        "risk.categ", "Freq", fill = "Histology", color="Histology", 
                        ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                        lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                               legend.title = element_text(size = 5), legend.text = element_text(size = 5))

ggsave(ggarrange(age.p, grade.p, X1p19q.p, IDH.p, histype.p, ncol=4, nrow=2), width = 20, height = 12, units = 'cm', 
       file='/result/section5/lgglike/plgg_325_cli_risk_com.pdf')



####################CGGA 301
sub.cli.data <- readRDS(file='/data/LggRiskScores/plgg_301_risk_score.rds')

# > mean(sub.cli.data$Age, na.rm = T)
# [1] 39.7451

sub.cli.data <- sub.cli.data %>% mutate(os = Censor..alive.0..dead.1., os_time = OS, 
                                        age_categ = ifelse(Age >= 40, '>= 40', '<40'), gender = factor(Gender,levels = c('Female', 'Male')), 
                                        TCGA_subtypes = factor(TCGA_subtypes, levels = c('Mesenchymal', 'Neural', 'Proneural', 'Classical')), 
                                        Histology = factor(Histology, levels = c('O', 'OA', 'AO', 'AOA', 'A', 'AA')), 
                                        Grade = factor(Grade, levels = c('WHO II', 'WHO III')),
                                        idh.status = factor(IDH_mutation_status, levels = c('Mutant', 'Wildtype')), 
                                        X1p19q.codeletion = factor(X1p19q_Codeletion_status, levels = c('Codel', 'Non-codel')), 
                                        mgmt.status = factor(MGMTp_methylation_status, levels = c('un-methylated', 'methylated')), 
                                        TMZ.treated = factor(Chemo_status..TMZ.treated.1.un.treated.0., levels = c(0, 1)), 
                                        Radio.status = factor(Radio_status..treated.1.un.treated.0., levels = c(0, 1)), 
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, age_categ, gender, TCGA_subtypes, Grade, idh.status, 
                                        X1p19q.codeletion, mgmt.status, Histology, risk.categ) # TMZ.treated, Radio.status



# Fisher's Exact Test
p.values <- sapply(colnames(cli.sig.char)[2:9], function(char){
  # print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  print(table(cli.sig.char[, c(char, 'risk.categ')]))
  # 
  # or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  stat <- table(cli.sig.char[, c(char, 'risk.categ')])
  p.value <- chisq.test(stat[rowSums(stat) > 0, ])$p.value
  # return(c(or, p.value))
  return(p.value)
})

names(p.values) <- colnames(cli.sig.char)[2:9]

# age_categ            gender     TCGA_subtypes             Grade        idh.status X1p19q.codeletion       mgmt.status         Histology 
# 1.672779e-02      8.945838e-01      2.924751e-10      2.705489e-06      8.934893e-03      6.293555e-02      6.226581e-01      8.907676e-05 

####Plot
library('ggpubr')
age.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('age_categ', 'risk.categ')])), 
                       "risk.categ", "Freq", fill = "age_categ", color="age_categ", 
                       ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                       lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                              legend.title = element_text(size = 5), legend.text = element_text(size = 5))

subtype.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('TCGA_subtypes', 'risk.categ')])), 
                   "risk.categ", "Freq", fill = "TCGA_subtypes", color="TCGA_subtypes", 
                   ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                   lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                          legend.title = element_text(size = 5), legend.text = element_text(size = 5))

grade.p <- ggbarplot(as.data.frame(table(cli.sig.char[, c('Grade', 'risk.categ')])), 
                     "risk.categ", "Freq", fill = "Grade", color="Grade", 
                     ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                     lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                            legend.title = element_text(size = 5), legend.text = element_text(size = 5))

IDH.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('idh.status', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "idh.status", color="idh.status", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))

histype.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('Histology', 'risk.categ')])), 
                        "risk.categ", "Freq", fill = "Histology", color="Histology", 
                        ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                        lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                               legend.title = element_text(size = 5), legend.text = element_text(size = 5))

ggsave(ggarrange(age.p, subtype.p, grade.p, IDH.p, histype.p, ncol=4, nrow=2), width = 20, height = 12, units = 'cm', 
       file='/result/section5/lgglike/plgg_301_cli_risk_com.pdf')




################################GSE107850
sub.cli.data <- readRDS(file='/data/LggRiskScores/plgg_GSE107850_risk_score.rds')
# > mean(sub.cli.data$age)
# [1] 43.94371
sub.cli.data <- sub.cli.data %>% mutate(age_categ = ifelse(age >= 44, '>= 44', '<44'), gender = factor(gender,levels = c('Female', 'Male')), 
                                        histology = factor(histology, levels = c('AOD GrII', 'AOA GrII', 'AA GrII')), 
                                        idh.status = factor(idh.status, levels = c('mutated', 'normal', 'undetermined')),
                                        performance = factor(performance, levels = c('0', '1', '2')),
                                        therapy = factor(therapy, levels = c('RT', 'TMZ')), 
                                        type.of.sugery = factor(type.of.sugery, levels = c('Total removal', 'Biopsy', 'Partial removal')), 
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, age_categ, gender, histology, idh.status, performance, therapy, type.of.sugery, risk.categ)



# Fisher's Exact Test
p.values <- sapply(colnames(cli.sig.char)[2:8], function(char){
  # print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  print(table(cli.sig.char[, c(char, 'risk.categ')]))
  # 
  # or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  stat <- table(cli.sig.char[, c(char, 'risk.categ')])
  p.value <- chisq.test(stat[rowSums(stat) > 0, ])$p.value
  # return(c(or, p.value))
  return(p.value)
})

names(p.values) <- colnames(cli.sig.char)[2:8]


# age_categ         gender      histology     idh.status    performance        therapy type.of.sugery 
# 3.025184e-01   1.057126e-01   2.178359e-02   6.077613e-07   2.260280e-01   1.000000e+00   1.676312e-01 
####Plot
library('ggpubr')
histology.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('histology', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "histology", color="histology", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))

IDH.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('idh.status', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "idh.status", color="idh.status", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))


ggsave(ggarrange(histology.p, IDH.p, ncol=4, nrow=1), width = 20, height = 6, units = 'cm', 
       file='/result/section5/lgglike/plgg_GSE107850_cli_risk_com.pdf')




################################GSE55918
sub.cli.data <- readRDS(file='/data/LggRiskScores/plgg_GSE55918_risk_score.rds')
# > mean(sub.cli.data$Age)
# [1] 44.26177
sub.cli.data <- sub.cli.data %>% mutate(age_categ = ifelse(Age >= 44, '>= 44', '<44'),
                                        os = Survival.time.months, os.event = Censored, 
                                        histology = factor(Histological, levels = c('Oligodendroglioma', 'Anaplastic oligodendroglioma', 
                                                                                    'Anaplastic astrocytoma', 'Astrocytoma', 'Mixed')), 
                                        grade = factor(Tumor.grade, levels = c('G2', 'G3')),
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, age_categ, histology, grade, risk.categ)



# Fisher's Exact Test
p.values <- sapply(colnames(cli.sig.char)[2:4], function(char){
  # print(rowSums(table(cli.sig.char[, c(char, 'risk.categ')])))
  print(table(cli.sig.char[, c(char, 'risk.categ')]))
  # 
  # or <- fisher.test(table(cli.sig.char[, c(char, 'risk.categ')]))$estimate
  stat <- table(cli.sig.char[, c(char, 'risk.categ')])
  p.value <- chisq.test(stat[rowSums(stat) > 0, ])$p.value
  # return(c(or, p.value))
  return(p.value)
})

names(p.values) <- colnames(cli.sig.char)[2:4]

# age_categ    histology        grade 
# 0.0005842569 0.0364849257 0.0044904304 

####Plot
library('ggpubr')
age.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('age_categ', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "age_categ", color="age_categ", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))

histology.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('histology', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "histology", color="histology", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
legend.title = element_text(size = 5), legend.text = element_text(size = 5))

grade.p  <- ggbarplot(as.data.frame(table(cli.sig.char[, c('grade', 'risk.categ')])), 
                    "risk.categ", "Freq", fill = "grade", color="grade", 
                    ylab='Number of patients', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
                    lab.vjust=2, palette="Paired") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                                                           legend.title = element_text(size = 5), legend.text = element_text(size = 5))

ggsave(ggarrange(age.p, grade.p, ncol=4, nrow=1), width = 20, height = 6, units = 'cm',  
       file='/result/section5/lgglike/plgg_GSE55918_cli_risk_com.pdf')




