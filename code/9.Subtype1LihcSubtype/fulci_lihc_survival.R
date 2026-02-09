
# load data
fulci.linc.cli.data <- readRDS(file='F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/data/lci_xinweiwang_cli_data.rds')
load(file='F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/data/lci_xinweiwang_expr_data.RData')
load(file='D:/CancerNeuroscience/Github/data/lihc_lasso_binomial_res.RData')


load('F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/data/exp_gene_anno.RData')
fulci.lihc.exp.data <- gene.max.exp.profile
fulci.lihc.exp.data <- fulci.lihc.exp.data %>% dplyr::select(-ID, -ENTREZ_GENE_ID, -`Gene Symbol`)
fulci.lihc.exp.data <- fulci.lihc.exp.data[intersect(rownames(fulci.lihc.exp.data), exp.gene.anno$GeneID), ]
rownames(fulci.lihc.exp.data) <- exp.gene.anno$Symbol[match(rownames(fulci.lihc.exp.data), exp.gene.anno$GeneID)]


# data prepare
fulci.linc.cli.data <- subset(fulci.linc.cli.data, Tissue.Type == 'Tumor')
fulci.linc.cli.data <- subset(fulci.linc.cli.data, !is.na(Survival.status))

# transform labels from months to days
fulci.linc.cli.data$Survival.months <- fulci.linc.cli.data$Survival.months*30.4375
fulci.linc.cli.data$Recurr.months <- fulci.linc.cli.data$Recurr.months*30.4375

fulci.lihc.exp.data <- fulci.lihc.exp.data[, fulci.linc.cli.data$Affy_GSM]


# risk evaluation
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
  
  # exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
  exp.dat <- as.matrix(exp.dat[gene.set, ])
  
  risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
  if(!is.null(cut.off)){
    risk.categ <- ifelse(risk.score >= cut.off, 'high risk', 'low risk')
    
  }else{
    risk.categ <- ifelse(risk.score >= quantile(risk.score, 0.60), 'high risk', 'low risk')
    
  }
  return(data.frame(risk.score, risk.categ))
}

active.k.vals <- subset(active.k.vals, symbol %in% rownames(fulci.lihc.exp.data)) # 39/43
lci.risk.score <- RiskEsti(exp.dat=fulci.lihc.exp.data, gene.set=active.k.vals$symbol, risk.coef=active.k.vals$coef, cut.off=NULL)
fulci.linc.cli.data <- cbind(fulci.linc.cli.data, lci.risk.score[fulci.linc.cli.data$Affy_GSM, ])



# survival plot 
source('F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/code/Rscript/survival_plot.R')

SurvivalPlot(survival.data=fulci.linc.cli.data[, c('LCS.ID', 'Survival.months', 'Survival.status')], 
             sample.class=fulci.linc.cli.data[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_os.pdf', 
             out.file.path='D:/CancerNeuroscience/Github/result/section5/lihclike/')

SurvivalPlot(survival.data=fulci.linc.cli.data[, c('LCS.ID', 'Recurr.months', 'Recurr.status')], 
             sample.class=fulci.linc.cli.data[, c('LCS.ID', 'risk.categ')], filename='fulci_lihc_rfs.pdf', 
             out.file.path='D:/CancerNeuroscience/Github/result/section5/lihclike/')


saveRDS(fulci.linc.cli.data, file='D:/CancerNeuroscience/Github/data/lihcRiskScores/fulci_risk_score.rds')


######COX
# load data
lci.risk.score <- readRDS(file='D:/CancerNeuroscience/Github/data/lihcRiskScores/fulci_risk_score.rds')

lci.risk.score <- lci.risk.score[, c(3, 5, 	8:18, 24, 19, 20)]
colnames(lci.risk.score) <- c('Affy_GSM',  'Metastasis_risk', 'Gender', 'Age', 'HBV_status', 'ALT', 'Tumor_size', 
                              'Multinodular', 'Cirrhosis', 'TNM_staging', 'BCLC_staging', 'CLIP_staging', 'AFP', 'risk.categ', 'os', 'os_time')

# data prepare
lci.risk.score$Metastasis_risk <- factor(lci.risk.score$Metastasis_risk, levels=c('low', 'high'))
lci.risk.score$Age <- ifelse(lci.risk.score$Age < 60, '<60', '≥60')
lci.risk.score$ALT <- factor(lci.risk.score$ALT, levels=c('low', 'high'))

lci.risk.score$Tumor_size[lci.risk.score$Tumor_size == '.'] <- NA
lci.risk.score$Tumor_size <- factor(lci.risk.score$Tumor_size, levels=c('small', 'large'))

lci.risk.score$TNM_staging[lci.risk.score$TNM_staging == '.'] <- NA
lci.risk.score$TNM_staging <- ifelse(lci.risk.score$TNM_staging %in% c('I', 'II'), 'I/II', 'III/IV')
lci.risk.score$BCLC_staging[lci.risk.score$BCLC_staging == '.'] <- NA
lci.risk.score$BCLC_staging <- ifelse(lci.risk.score$BCLC_staging %in% c('0', 'A'), '0/A', 'B/C')
lci.risk.score$CLIP_staging[lci.risk.score$CLIP_staging == '.'] <- NA
lci.risk.score$CLIP_staging <- ifelse(lci.risk.score$CLIP_staging %in% c('0', '1'), '0/1', '≥2')
lci.risk.score$AFP[lci.risk.score$AFP == '.'] <- NA

lci.risk.score$AFP <- factor(lci.risk.score$AFP, levels=c('low', 'high'))
lci.risk.score$risk.categ <- factor(lci.risk.score$risk.categ, levels=c('low risk', 'high risk'))


cli.sig.char <- lci.risk.score


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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[2:14])
#            variate                     HR p.value
# risk.categhigh risk 1.6347 (1.0676-2.5031)  0.0238

# 多因素cox分析
cli.sig.char <- cli.sig.char[, c('Affy_GSM', 'os', 'os_time', 'risk.categ', 'Age', 'Gender', 
                                 'Cirrhosis', 'AFP', "TNM_staging", "BCLC_staging", "CLIP_staging")]
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]

source('F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/code/Rscript/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:8, 11))

# multiv HR (95% CI for HR) multiv p value
# 0.0241    1.7164 (1.1097-2.6547)         0.0152

















