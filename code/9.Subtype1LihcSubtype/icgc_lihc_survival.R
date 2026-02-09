# load data
load(file='/data/lihc_lasso_binomial_res.RData')
icgc.linc.cli.data <- readRDS(file='/data/icgc_linc_cli_data.rds')
icgc.lihc.vst.exp <- readRDS(file='/data/icgc_vst_norm_geneExp.rds')


load('/data/exp_gene_anno.RData')
icgc.lihc.vst.exp <- icgc.lihc.vst.exp[as.character(subset(exp.gene.anno, !is.na(ICGCExpGeneFilter))$GeneID), ]
rownames(icgc.lihc.vst.exp) <- exp.gene.anno$Symbol[match(rownames(icgc.lihc.vst.exp), exp.gene.anno$GeneID)]


# data prepare
# active.k.vals <- active.k.vals.1se
icgc.linc.cli.data <- subset(icgc.linc.cli.data, !filter)
rownames(icgc.linc.cli.data) <- icgc.linc.cli.data$icgc_sample_id


# risk evaluation
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
  
  # exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
  exp.dat <- as.matrix(exp.dat[gene.set, ])
  
  risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
  if(!is.null(cut.off)){
    risk.categ <- ifelse(risk.score >= cut.off, 'high risk', 'low risk')
    
  }else{
    risk.categ <- ifelse(risk.score >= quantile(risk.score, 0.6), 'high risk', 'low risk')
    
  }
  return(data.frame(risk.score, risk.categ))
}


# all patients
active.k.vals <- subset(active.k.vals, symbol %in% rownames(icgc.lihc.vst.exp)) # 40/43
icgc.risk.score <- RiskEsti(exp.dat=icgc.lihc.vst.exp, gene.set=active.k.vals$symbol, risk.coef=active.k.vals$coef, cut.off=NULL)
icgc.linc.cli.data <- merge(icgc.linc.cli.data, icgc.risk.score, by='row.names')



# survival plot 
source('/code/Rscript/survival_plot.R')
icgc.linc.cli.data$donor_vital_status <- ifelse(icgc.linc.cli.data$donor_vital_status == 'deceased', 1, 0)
SurvivalPlot(survival.data=icgc.linc.cli.data[, c('icgc_donor_id', 'donor_survival_time', 'donor_vital_status')], 
             sample.class=icgc.linc.cli.data[, c('icgc_donor_id', 'risk.categ')], filename='icgc_lihc_os.pdf', 
             out.file.path='/result/section5/lihclike/')

saveRDS(icgc.linc.cli.data, file='/data/lihcRiskScores/icgc_risk_score.rds')

##########Cox
# load data
icgc.lihc.risk.score <- readRDS(file='/data/lihcRiskScores/icgc_risk_score.rds')
icgc.lihc.cli.data <- readRDS(file='/data/icgc_linc_cli_data.rds')

# data prepare
sig.score <- icgc.lihc.risk.score[, c('icgc_donor_id', 'donor_vital_status', 'donor_survival_time', 'risk.score', 'risk.categ')]
colnames(sig.score) <- c('icgc_donor_id', 'os', 'os_time', 'risk.score', 'risk.categ')
sig.score$risk.categ <- factor(sig.score$risk.categ, levels=c('low risk', 'high risk'))


icgc.lihc.cli.data$tumor_size_category <- ifelse(icgc.lihc.cli.data$Tumor.size..mm. > 50, '>50', '≤50')
icgc.lihc.cli.data$tumor_size_category <- factor(icgc.lihc.cli.data$tumor_size_category, levels=c('≤50', '>50'))

cli.char <- icgc.lihc.cli.data[, c('icgc_donor_id', 'age_category', 'gender_category', 'Smoking', 
                                   'tumor_size_category', 'virus_infection_category', 'pathologic_stage_category', 'histologic_grade_category', 
                                   'portal_invasion_category', 'hepatic_invasion_category', 'bile_invasion_category', 'fibrosisc_category')]

cli.sig.char <- merge(sig.score, cli.char, by='icgc_donor_id')

# 单因素cox分析
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
    p.value <- round(tmp$coefficients[ ,5], digits = 4);
    p.value[which(p.value < 0.0001)] <- "<0.0001";
    
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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[5:16])
#             variate                     HR p.value
# risk.categhigh risk 1.9897 (1.0223-3.8726)  0.0429

# 多因素cox分析
cli.sig.char <- cli.sig.char[, c('icgc_donor_id', 'os', 'os_time', 'risk.categ', 'age_category', 'gender_category', 
                                 'pathologic_stage_category', 'histologic_grade_category', 'fibrosisc_category')]
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]


source('/code/Rscript/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:9))

# multiv HR (95% CI for HR) multiv p value
# 1.522 (0.7564-3.0626)         0.2391


