
# load data
load(file='D:/CancerNeuroscience/Github/data/lihc_lasso_binomial_res.RData')
load(file='F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/data/lec_snorri_s_thorgeirsson_expr_data.RData')
nci.linc.cli.data <- readRDS(file='F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/data/lec_snorri_s_thorgeirsson_cli_data.rds')


load('F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/data/exp_gene_anno.RData')
nci.lihc.exp.data <- batch.expr.max.data
nci.lihc.exp.data <- nci.lihc.exp.data %>% dplyr::select(-GeneID, -Symbol, -ID)
nci.lihc.exp.data <- nci.lihc.exp.data[intersect(rownames(nci.lihc.exp.data), exp.gene.anno$GeneID), ]
rownames(nci.lihc.exp.data) <- exp.gene.anno$Symbol[match(rownames(nci.lihc.exp.data), exp.gene.anno$GeneID)]


# data prepare
nci.linc.cli.data <- subset(nci.linc.cli.data, !is.na(OS_Status))
nci.lihc.exp.data <- nci.lihc.exp.data[, intersect(nci.linc.cli.data$Array, colnames(nci.lihc.exp.data))]
nci.linc.cli.data <- subset(nci.linc.cli.data, Array %in% colnames(nci.lihc.exp.data))


# transform labels from months to days
nci.linc.cli.data$OS_Time <- nci.linc.cli.data$OS_Time*30.4375
rownames(nci.linc.cli.data) <- nci.linc.cli.data$Array

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


active.k.vals <- subset(active.k.vals, symbol %in% rownames(nci.lihc.exp.data)) # 41/43
lec.risk.score <- RiskEsti(exp.dat=nci.lihc.exp.data, gene.set=active.k.vals$symbol, risk.coef=active.k.vals$coef, cut.off=NULL)
nci.linc.cli.data <- cbind(nci.linc.cli.data, lec.risk.score[nci.linc.cli.data$Array, ])


# survival plot 
source('F:/PostdoctoralDataBackup/DesktopCP/E/ExpressionITHProject/GitHub/code/Rscript/survival_plot.R')

nci.linc.cli.data$OS_Status[nci.linc.cli.data$OS_Time > 365.25 *5] <- 0
nci.linc.cli.data$OS_Time[nci.linc.cli.data$OS_Time > 365.25 *5] <- 365.25 * 5

SurvivalPlot(survival.data=nci.linc.cli.data[, c('Array', 'OS_Time', 'OS_Status')], 
             sample.class=nci.linc.cli.data[, c('Array', 'risk.categ')], filename='nci_lihc_os.pdf', 
             out.file.path='D:/CancerNeuroscience/Github/result/section5/lihclike/')

saveRDS(nci.linc.cli.data, file='D:/CancerNeuroscience/Github/data/lihcRiskScores/nci_risk_score.rds')

########### Cox
cli.sig.char <- readRDS(file='D:/CancerNeuroscience/Github/data/lihcRiskScores/nci_risk_score.rds')
colnames(cli.sig.char)[2:3] <- c('os_time','os')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))


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
    p.value <- tmp$coefficients[ ,5];
    # p.value[which(p.value < 0.0001)] <- "<0.0001";
    
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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[6])
#                       variate                     HR            p.value
# risk.categ risk.categhigh risk 1.7891 (1.1029-2.9022) 0.0184259918429043
