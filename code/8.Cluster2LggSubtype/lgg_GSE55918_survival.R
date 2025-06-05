
# clinical data
cli.data <- read.csv(file = '/data/GliomaData/GSE55918_clinical_data.txt', header = T, sep = '\t', comment.char = '#')


cli.data <- cli.data %>% dplyr::rename(Tumor.grade = characteristics..Tumor.grade.diagnosis, Histological = characteristics..Histological.dianosis,
                                Survival.time.months = characteristics..Survival.time..months., 
                                Age = characteristics..Age.at.diag1sis..years., Censored = characteristics..Censored, 
                                Treatment.type = characteristics..Treatment.type, Chemotherapy.drug.name = characteristics..Chemotherapy.drug.name)



cli.data <- cli.data %>% mutate(Tumor.grade = trimws(Tumor.grade), Histological = trimws(Histological), 
                                Survival.time.months = ifelse(Survival.time.months == 'Null', NA, Survival.time.months),
                                Age = ifelse(Age == 'Null', NA, Age), Treatment.type = ifelse(Treatment.type == 'Null', NA, Treatment.type),
                                Chemotherapy.drug.name = ifelse(Chemotherapy.drug.name == 'Null', NA, Chemotherapy.drug.name)) %>% 
  mutate(Age = as.numeric(Age), Survival.time.months = as.numeric(Survival.time.months), 
         Tumor.grade = ifelse(Tumor.grade == 'Null', NA, Tumor.grade), Histological = ifelse(Histological == 'Null', NA, Histological))


cli.data <- cli.data %>% filter(Tumor.grade %in% c('G2', 'G3'), data.set.name !='TCGA', !is.na(Survival.time.months))
cli.data <- cli.data %>% remove_rownames() %>% column_to_rownames(var = 'raw.data.file')

# exp
norm.exp <- read.csv(file = '/data/GliomaData/GSE55918_Matrix_GliomaClusteringAnalysis.txt', header = T, sep = '\t')
norm.exp <- norm.exp %>% column_to_rownames(var = 'Probe.set.ID')
norm.exp <- norm.exp[, cli.data$raw.data.file]


probe.anno1 <- readxl::read_xls("/data/GliomaData/BA1_probe annotation.xls", col_names = T)
probe.anno2 <- readxl::read_xls("/data/GliomaData/BA2_probe annotation.xls", col_names = T)
colnames(probe.anno2)[2] <- 'Symbol'
probe.anno1 <- probe.anno1 %>% column_to_rownames(var = 'orginal_id')
probe.anno2 <- probe.anno2 %>% column_to_rownames(var = 'probe')


norm.exp1 <- merge(probe.anno1[, 'symbols', FALSE], norm.exp, by = 'row.names')
norm.exp2 <- merge(probe.anno2, norm.exp, by = 'row.names')

norm.exp1 <- subset(norm.exp1, symbols %in% active.k.vals.1se$symbol)
norm.exp1 <- norm.exp1[, -c(1)] %>% remove_rownames() %>% column_to_rownames(var = 'symbols')
norm.exp2 <- subset(norm.exp2, Symbol %in% active.k.vals.1se$symbol)
norm.exp2 <-norm.exp2[, -c(1)] %>% remove_rownames() %>% column_to_rownames(var = 'Symbol')

ana.data <- rbind(norm.exp1, norm.exp2)
ana.data <- distinct(ana.data)




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


risk.score <- RiskEsti(exp.dat=ana.data, gene.set=active.k.vals.1se$symbol, risk.coef=active.k.vals.1se$coef, cut.off='Top')
sub.cli.data <- merge(cli.data, risk.score, by='row.names')
colnames(sub.cli.data)[1] <- 'patient_id'


# survival plot 
source('/code/survival_plot.R')

sub.cli.data$Survival.time.months <- round(sub.cli.data$Survival.time.months*30.25)
SurvivalPlot(survival.data=sub.cli.data[, c('patient_id', 'Survival.time.months', 'Censored')], 
             sample.class=sub.cli.data[, c('patient_id', 'risk.categ')], filename='plgg_GSE55918_os.pdf', 
             out.file.path='/result/section6/lgglike/')


saveRDS(sub.cli.data, file='/data/LggRiskScores/plgg_GSE55918_risk_score.rds')




################################uni_mul_cox
sub.cli.data <- readRDS(file='/data/LggRiskScores/plgg_GSE55918_risk_score.rds')
# > mean(sub.cli.data$Age)
# [1] 44.26177
sub.cli.data <- sub.cli.data %>% mutate(age_categ = ifelse(Age >= 44, '>= 44', '<44'),
                                        os = Survival.time.months, os.event = Censored, 
                                        histology = factor(Histological, levels = c('Oligodendroglioma', 'Anaplastic oligodendroglioma', 
                                                                                    'Anaplastic astrocytoma', 'Astrocytoma', 'Mixed')), 
                                        grade = factor(Tumor.grade, levels = c('G2', 'G3')),
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% dplyr::select(patient_id, os, os.event, age_categ, histology, grade, risk.categ)


UnivariateCox <- function(cli.data, covariates)
{
  library('survival')
  #STEP1:构建单因素分析的对象
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(os, os.event)~', x)));
  
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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[4:7])




# 多因素cox分析
cli.sig.char <- subset(cli.sig.char, idh.status != 'undetermined')
# source('/code/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os, event=cli.sig.char$os.event, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:7))


