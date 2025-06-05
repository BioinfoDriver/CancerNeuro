library('dplyr')
###################
mc3MutData <- readRDS(file = '/data/mc3MutData.rds')

mc3MutData <- subset(mc3MutData, Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',  'In_Frame_Ins', 
                                                               'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'))


cbiPortalData <- subset(mc3MutData, Hugo_Symbol %in% c('CHRM3', 'HTR5A', 'GABRG1'), 
                        select = c(SAMPLE_BARCODE, DISEASE, Hugo_Symbol, HGVSp_Short, Variant_Classification, 
                                   Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, t_alt_count, t_ref_count))
cbiPortalData <- cbiPortalData %>% mutate(HGVSp_Short = gsub('p.', '', HGVSp_Short))

colnames(cbiPortalData) <- c('Sample_ID',	'Cancer_Type', 'Hugo_Symbol', 'Protein_Change', 'Mutation_Type', 'Chromosome',	
                             'Start_Position', 'End_Position', 'Reference_Allele', 'Variant_Allele', 'Tumor_Alt_Count', 'Tumor_Ref_Count')


write.table(cbiPortalData, col.names = T, row.names = F, sep = '\t', quote = F, file = '/result/section1/cbiPortalData.txt')
