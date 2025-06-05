
library('dplyr')
# read count data
norm.exp.data <- read.csv(file = '/data/GliomaData/CGGA.mRNA_array_301_gene_level.20200506.txt', 
                      header = T, sep = '\t', stringsAsFactors = F)

norm.exp.data <- norm.exp.data %>% tibble::column_to_rownames(var = 'Gene_Name')

# read clinical data
cli.data <- read.csv(file = '/data/GliomaData/CGGA.mRNA_array_301_clinical.20200506.txt', 
                     header = T, sep = '\t', stringsAsFactors = F)

cli.data <- cli.data %>% tibble::column_to_rownames(var = 'CGGA_ID')
cli.data <- cli.data[colnames(norm.exp.data), ]

cli.data <- cli.data %>% mutate(disType = case_when(PRS_type == 'Primary' & Histology == 'GBM' ~ 'pGBM', 
                                                    PRS_type == 'Recurrent' & Histology == 'rGBM' ~ 'rGBM',
                                                    PRS_type == 'Primary' & Histology %in% c('A', 'AA', 'AO', 'AOA', 'O', 'OA') ~ 'pLGG',
                                                    PRS_type == 'Recurrent' & Histology %in% c('rA', 'rAA', 'rAO', 'rAOA', 'rO', 'rOA') ~ 'rLGG',
                                                    PRS_type == 'Secondary' & Histology == 'sGBM' ~ 'sGBM'))


sub.cli.data <- subset(cli.data, disType == 'pLGG')

exp.data <- norm.exp.data[, intersect(colnames(norm.exp.data), rownames(sub.cli.data))]
sub.cli.data <- sub.cli.data[intersect(colnames(exp.data), rownames(sub.cli.data)), ]



# load data
load(file='/data/lgg_lasso_binomial_res.RData')


# risk evaluation
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
  
  # exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
  exp.dat <- as.matrix(exp.dat[gene.set, ])
  
  risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
  if(!is.null(cut.off)){
    risk.categ <- ifelse(risk.score >= quantile(risk.score, 0.8), 'high risk', 'low risk')
    
  }else{
    risk.categ <- ifelse(risk.score >= median(risk.score), 'high risk', 'low risk')
    
  }
  return(data.frame(risk.score, risk.categ))
}


risk.score <- RiskEsti(exp.dat=exp.data, gene.set=active.k.vals.1se$symbol, risk.coef=active.k.vals.1se$coef, cut.off='Top')
sub.cli.data <- merge(sub.cli.data, risk.score, by='row.names')
colnames(sub.cli.data)[1] <- 'patient_id'


# survival plot 
source('/code/survival_plot.R')


SurvivalPlot(survival.data=sub.cli.data[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=sub.cli.data[, c('patient_id', 'risk.categ')], filename='plgg_301_os.pdf', 
             out.file.path='/result/section6/lgglike/')


saveRDS(sub.cli.data, file='/data/LggRiskScores/plgg_301_risk_score.rds')



################################uni_mul_cox
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

cli.sig.char <- sub.cli.data %>% select(patient_id, os_time, os, age_categ, gender, TCGA_subtypes, Grade, idh.status, 
                                        risk.categ, mgmt.status, Histology) # X1p19q.codeletion, TMZ.treated, Radio.status

cli.sig.char <- subset(cli.sig.char, !is.na(os_time) & !is.na(os))

UnivariateCox <- function(cli.data, covariates)
{
  library('survival')
  #STEP1:构建单因素分析的对象
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(os_time, os)~', x)));
  
  #STEP2:单因素Cox分析
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = cli.data)});
  
  #STEP3:提取有用信息
  univ_results <- lapply(univ_models, function(x)
  {                             
    tmp <-summary(x);
    
    #提取p值，保留两位有效数字
    # p.value <- round(tmp$coefficients[ ,5], digits = 4);
    # p.value[which(p.value < 0.0001)] <- "<0.0001";
    p.value <- tmp$coefficients[ ,5]
    
    #提取beta值，这里的coefficients为矩阵，但是只有一行，所以可以这样取值
    #beta <- round(tmp$coefficients[ ,1], digits = 4);
    
    #提取风险比
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    
    #提取95%置信区间上下界
    HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 4);
    HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 4);    
    
    #合并风险比HR和置信区间为一个内容
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
    
    variate <- rownames(tmp$coefficients);
    
    #将所有值合并在一个矩阵中
    all.data <- as.data.frame(cbind(variate, HR, p.value));
  }
  )
  univ_results <- do.call(rbind, univ_results)
  return(univ_results)
}

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[4:11])




# 多因素cox分析
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]
# cli.sig.char <- subset(cli.sig.char, idh.status != 'undetermined')
# source('/code/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:11))
