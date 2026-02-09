

nReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
load(file='/data/lgg_lasso_binomial_res.RData')
S4signature <- active.k.vals.1se$entrez_id
S4signature <- subset(nReceptors, NCBI.Gene.ID %in% S4signature)$Approved.symbol

NTRLGGindex <- c('HTR5A','HTR3B','HTR2A','HTR1E','GRM7','GRM4','GRM3','GRM2','GRIN2B',
'GRIN2A','GRIN1','GRIK4','GRID2','GABRG2','GABRB2','GABRA5','GABRA4',
'GABRA2','GABRA1','DRD1','CHRNB1','CHRNA2','CHRM3','CHRM1','ADRA1B')


> intersect(NTRLGGindex, S4signature)
[1] "GRIN2B" "CHRNB1"


# table(subset(nReceptors, Approved.symbol %in% S4signature)$classesOfNeuroreceptors)
# CholinergicReceptor        GABAReceptor   GlutamateReceptor     GlycineReceptor   HistamineReceptor  PurinergicReceptor 
#                   4                   3                   7                   1                   1                   2

table(subset(nReceptors, Approved.symbol %in% NTRLGGindex)$classesOfNeuroreceptors)
# AdrenergicReceptor CholinergicReceptor    DopamineReceptor        GABAReceptor   GlutamateReceptor   SerotoninReceptor 
#                  1                   4                   1                   6                   9                   4 

# > subset(nReceptors, Approved.symbol %in% S4signature) %>% select(Approved.symbol, classesOfNeuroreceptors) %>% arrange(desc(classesOfNeuroreceptors))
# Approved.symbol classesOfNeuroreceptors
# ADORA2A         ADORA2A      PurinergicReceptor
# P2RY2             P2RY2      PurinergicReceptor
# HRH1               HRH1       HistamineReceptor
# GLRA3             GLRA3         GlycineReceptor
# GRIA2             GRIA2       GlutamateReceptor
# GRID1             GRID1       GlutamateReceptor
# GRIK1             GRIK1       GlutamateReceptor
# GRIN2B           GRIN2B       GlutamateReceptor
# GRIN3A           GRIN3A       GlutamateReceptor
# GRM1               GRM1       GlutamateReceptor
# GRM6               GRM6       GlutamateReceptor
# GABRG1           GABRG1            GABAReceptor
# GABRQ             GABRQ            GABAReceptor
# GABBR1           GABBR1            GABAReceptor
# CHRNA1           CHRNA1     CholinergicReceptor
# CHRNA4           CHRNA4     CholinergicReceptor
# CHRNA5           CHRNA5     CholinergicReceptor
# CHRNB1           CHRNB1     CholinergicReceptor

# > subset(nReceptors, Approved.symbol %in% NTRLGGindex) %>% select(Approved.symbol, classesOfNeuroreceptors) %>% arrange(desc(classesOfNeuroreceptors))
# Approved.symbol classesOfNeuroreceptors
# HTR1E            HTR1E       SerotoninReceptor
# HTR2A            HTR2A       SerotoninReceptor
# HTR5A            HTR5A       SerotoninReceptor
# HTR3B            HTR3B       SerotoninReceptor
# GRID2            GRID2       GlutamateReceptor
# GRIK4            GRIK4       GlutamateReceptor
# GRIN1            GRIN1       GlutamateReceptor
# GRIN2A          GRIN2A       GlutamateReceptor
# GRIN2B          GRIN2B       GlutamateReceptor
# GRM2              GRM2       GlutamateReceptor
# GRM3              GRM3       GlutamateReceptor
# GRM4              GRM4       GlutamateReceptor
# GRM7              GRM7       GlutamateReceptor
# GABRA1          GABRA1            GABAReceptor
# GABRA2          GABRA2            GABAReceptor
# GABRA4          GABRA4            GABAReceptor
# GABRA5          GABRA5            GABAReceptor
# GABRB2          GABRB2            GABAReceptor
# GABRG2          GABRG2            GABAReceptor
# DRD1              DRD1        DopamineReceptor
# CHRM1            CHRM1     CholinergicReceptor
# CHRM3            CHRM3     CholinergicReceptor
# CHRNA2          CHRNA2     CholinergicReceptor
# CHRNB1          CHRNB1     CholinergicReceptor
# ADRA1B          ADRA1B      AdrenergicReceptor

library(ggvenn)
library(ggplot2)
vennPlot <- ggvenn(list('NTR-LGG index' = NTRLGGindex, 'S4-like LGG signature' = S4signature), 
                   show_elements = T, set_name_size = 5, text_size = 3)

ggsave(filename = '/result/section6/lgglike/S4signatureNTRLGGindex.pdf', 
       plot = vennPlot, units = 'cm', width = 8, height = 8)



