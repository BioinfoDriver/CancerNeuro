
library('dplyr')
library('GEOquery')


# read exp data
gse <- getGEO(filename = '/data/GliomaData/GSE107850_series_matrix.txt')
norm.exp.data <- exprs(gse)
cli.data <- pData(gse) 
cli.data <- cli.data[, 42:51]
colnames(cli.data) <- c('age', 'gender', 'histology', 'idh.status', 'igs', 'performance', 'pfs', 'pfs.event', 'therapy', 'type.of.sugery')
cli.data <- cli.data %>% mutate(age = as.numeric(age), pfs = as.numeric(pfs))



# annotation data
anno.data <- read.csv(file = '/data/GliomaData/GPL14951-11332.txt', 
                     header = T, sep = '\t', stringsAsFactors = F, comment.char = '#')

anno.data <- subset(anno.data, Source == 'RefSeq')
anno.data <- anno.data %>% select(ID, Entrez_Gene_ID) %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'ID')


# Select the probe with the highest normalized intensity averaged over all samples
norm.exp.data <- merge(anno.data, norm.exp.data, by='row.names')
norm.exp.data <- tibble::column_to_rownames(norm.exp.data, var = 'Row.names')

norm.exp.data <- split.data.frame(x = norm.exp.data, f = ~Entrez_Gene_ID)

# probe with maximal expression
gene.max.exp <- lapply(norm.exp.data, function(exp.data){
  index <- which.max(rowSums(exp.data[, -1]))
  exp.data <- exp.data[index, , FALSE]
})

gene.max.exp <- do.call(rbind, gene.max.exp)
gene.max.exp <- gene.max.exp  %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'Entrez_Gene_ID')


sub.cli.data <- cli.data
norm.exp.data <- gene.max.exp

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


risk.score <- RiskEsti(exp.dat=exp.data, gene.set=active.k.vals.1se$entrez_id, risk.coef=active.k.vals.1se$coef, cut.off='Top')
sub.cli.data <- merge(sub.cli.data, risk.score, by='row.names')
colnames(sub.cli.data)[1] <- 'patient_id'


# survival plot 
source('/code/survival_plot.R')

sub.cli.data <- sub.cli.data %>% mutate(pfs.event = ifelse(pfs.event == 'Yes', TRUE, FALSE))

SurvivalPlot(survival.data=sub.cli.data[, c('patient_id', 'pfs', 'pfs.event')], 
             sample.class=sub.cli.data[, c('patient_id', 'risk.categ')], filename='plgg_GSE107850_pfs.pdf', 
             out.file.path='/result/section6/lgglike/')


saveRDS(sub.cli.data, file='/data/LggRiskScores/plgg_GSE107850_risk_score.rds')



################################uni_mul_cox
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

cli.sig.char <- sub.cli.data %>% select(patient_id, pfs, pfs.event, age_categ, gender, histology, idh.status, performance, therapy, type.of.sugery, risk.categ)


UnivariateCox <- function(cli.data, covariates)
{
  library('survival')
  #STEP1:构建单因素分析的对象
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(pfs, pfs.event)~', x)));
  
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
cli.sig.char <- subset(cli.sig.char, idh.status != 'undetermined')
# source('/code/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$pfs, event=cli.sig.char$pfs.event, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:11))



