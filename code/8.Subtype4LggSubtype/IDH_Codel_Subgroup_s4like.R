
#########################TCGA
tcga.s4.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
tcga_glioma_cli_mol <- readRDS(file = 'D:/CancerNeuroscience/data/GliomaData/tcga_glioma_cli_mol.rds')
tcga_glioma_cli_mol <- tcga_glioma_cli_mol %>% mutate(bcr_patient_barcode = paste0(bcr_patient_barcode, '-01'))

tcga.s4.data <- merge(tcga.s4.data, tcga_glioma_cli_mol[, c('bcr_patient_barcode', 'IDH_STATUS', 
                                                            'MGMT_PROMOTER_STATUS', 'IDH_1P19Q_SUBTYPE', 'IDH_CODEL_SUBTYPE')], 
                      by.x = 'patient_id', by.y = 'bcr_patient_barcode')

# survival plot 
source('/code/0.DataPreparation/survival_plot.R')
SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_STATUS == 'WT')[, c('patient_id', 'os_time', 'os')], 
             sample.class=subset(tcga.s4.data, IDH_STATUS == 'WT')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_IDHwt_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_STATUS == 'Mutant')[, c('patient_id', 'os_time', 'os')], 
             sample.class=subset(tcga.s4.data, IDH_STATUS == 'Mutant')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_IDHmt_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Non-codel')[, c('patient_id', 'os_time', 'os')], 
             sample.class=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Non-codel')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_Noncodel_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Codel')[, c('patient_id', 'os_time', 'os')], 
             sample.class=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Codel')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_Codel_os.pdf', 
             out.file.path='/result/section5/lgglike/')


SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_CODEL_SUBTYPE == 'IDHmut-non-codel')[, c('patient_id', 'os_time', 'os')], 
             sample.class=subset(tcga.s4.data, IDH_CODEL_SUBTYPE == 'IDHmut-non-codel')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_IDHmutnoncodel_os.pdf', 
             out.file.path='/result/section5/lgglike/')



SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_STATUS == 'WT')[, c('patient_id', 'pfi_time', 'pfi')], 
             sample.class=subset(tcga.s4.data, IDH_STATUS == 'WT')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_IDHwt_pfi.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_STATUS == 'Mutant')[, c('patient_id', 'pfi_time', 'pfi')], 
             sample.class=subset(tcga.s4.data, IDH_STATUS == 'Mutant')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_IDHmt_pfi.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Non-codel')[, c('patient_id', 'pfi_time', 'pfi')], 
             sample.class=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Non-codel')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_Noncodel_pfi.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Codel')[, c('patient_id', 'pfi_time', 'pfi')], 
             sample.class=subset(tcga.s4.data, IDH_1P19Q_SUBTYPE == 'Codel')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_Codel_pfi.pdf', 
             out.file.path='/result/section5/lgglike/')


SurvivalPlot(survival.data=subset(tcga.s4.data, IDH_CODEL_SUBTYPE == 'IDHmut-non-codel')[, c('patient_id', 'pfi_time', 'pfi')], 
             sample.class=subset(tcga.s4.data, IDH_CODEL_SUBTYPE == 'IDHmut-non-codel')[, c('patient_id', 'risk.categ')], filename='tcga_lgg_IDHmutnoncodel_pfi.pdf', 
             out.file.path='/result/section5/lgglike/')


#########################CGGA-mRNAseq1
CGGAmRNAseq1.s4.data <- readRDS(file='/data/LggRiskScores/plgg_693_risk_score.rds')

SurvivalPlot(survival.data=subset(CGGAmRNAseq1.s4.data, IDH_mutation_status == 'Wildtype')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq1.s4.data, IDH_mutation_status == 'Wildtype')[, c('patient_id', 'risk.categ')], filename='plgg_693_IDHwt_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAmRNAseq1.s4.data, IDH_mutation_status == 'Mutant')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq1.s4.data, IDH_mutation_status == 'Mutant')[, c('patient_id', 'risk.categ')], filename='plgg_693_IDHmut_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAmRNAseq1.s4.data, X1p19q_codeletion_status == 'Codel')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq1.s4.data, X1p19q_codeletion_status == 'Codel')[, c('patient_id', 'risk.categ')], filename='plgg_693_Codel_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAmRNAseq1.s4.data, X1p19q_codeletion_status == 'Non-codel')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq1.s4.data, X1p19q_codeletion_status == 'Non-codel')[, c('patient_id', 'risk.categ')], filename='plgg_693_Noncodel_os.pdf', 
             out.file.path='/result/section5/lgglike/')


#########################CGGA-mRNAseq2
CGGAmRNAseq2.s4.data <- readRDS(file='/data/LggRiskScores/plgg_325_risk_score.rds')

SurvivalPlot(survival.data=subset(CGGAmRNAseq2.s4.data, IDH_mutation_status == 'Wildtype')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq2.s4.data, IDH_mutation_status == 'Wildtype')[, c('patient_id', 'risk.categ')], filename='plgg_325_IDHwt_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAmRNAseq2.s4.data, IDH_mutation_status == 'Mutant')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq2.s4.data, IDH_mutation_status == 'Mutant')[, c('patient_id', 'risk.categ')], filename='plgg_325_IDHmut_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAmRNAseq2.s4.data, X1p19q_codeletion_status == 'Codel')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq2.s4.data, X1p19q_codeletion_status == 'Codel')[, c('patient_id', 'risk.categ')], filename='plgg_325_Codel_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAmRNAseq1.s4.data, X1p19q_codeletion_status == 'Non-codel')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAmRNAseq1.s4.data, X1p19q_codeletion_status == 'Non-codel')[, c('patient_id', 'risk.categ')], filename='plgg_325_Noncodel_os.pdf', 
             out.file.path='/result/section5/lgglike/')


#########################CGGA-mRNA-array1
CGGAarray1.s4.data <- readRDS(file='/data/LggRiskScores/plgg_301_risk_score.rds')

SurvivalPlot(survival.data=subset(CGGAarray1.s4.data, IDH_mutation_status == 'Wildtype')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAarray1.s4.data, IDH_mutation_status == 'Wildtype')[, c('patient_id', 'risk.categ')], filename='plgg_301_IDHwt_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(CGGAarray1.s4.data, IDH_mutation_status == 'Mutant')[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=subset(CGGAarray1.s4.data, IDH_mutation_status == 'Mutant')[, c('patient_id', 'risk.categ')], filename='plgg_301_IDHmut_os.pdf', 
             out.file.path='/result/section5/lgglike/')


#########################GEO-mRNA-array3
GEOarray3.s4.data <- readRDS(file='/data/LggRiskScores/plgg_GSE107850_risk_score.rds')

SurvivalPlot(survival.data=subset(GEOarray3.s4.data, idh.status == 'normal')[, c('patient_id', 'pfs', 'pfs.event')], 
             sample.class=subset(GEOarray3.s4.data, idh.status == 'normal')[, c('patient_id', 'risk.categ')], filename='plgg_GSE107850_IDHwt_pfs.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=subset(GEOarray3.s4.data, idh.status == 'mutated')[, c('patient_id', 'pfs', 'pfs.event')], 
             sample.class=subset(GEOarray3.s4.data, idh.status == 'mutated')[, c('patient_id', 'risk.categ')], filename='plgg_GSE107850_IDHmut_pfs.pdf', 
             out.file.path='/result/section5/lgglike/')


