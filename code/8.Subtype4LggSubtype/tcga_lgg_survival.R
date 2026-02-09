library('ggpubr')
library('ggplot2')
library('patchwork')
library('pROC')

# load data
tcga.cli.data <- readRDS(file = '/data/tcgaPanCanCliData.rds')
tcga.lgg.cli.data <- subset(tcga.cli.data, cancer_type == 'LGG')
rownames(tcga.lgg.cli.data) <- paste0(rownames(tcga.lgg.cli.data), '-01')

load(file = '/data/panCanGeneExpData.RData')
tcgaExpData <- panCanTurGeneExp

tcga.lgg.exp <- tcgaExpData[, intersect(colnames(tcgaExpData), rownames(tcga.lgg.cli.data))]
tcga.lgg.cli.data <- tcga.lgg.cli.data[intersect(rownames(tcga.lgg.cli.data), colnames(tcga.lgg.exp)), ]

# risk evaluation
load(file='/data/lgg_lasso_binomial_res.RData')

RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
  
  # exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
  exp.dat <- as.matrix(exp.dat[gene.set, ])
  exp.dat[is.na(exp.dat)] <- 0
  
  risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
  if(!is.null(cut.off)){
    risk.categ <- ifelse(risk.score >= quantile(risk.score, 0.80), 'high risk', 'low risk')
    
  }else{
    risk.categ <- ifelse(risk.score >= median(risk.score), 'high risk', 'low risk')
    
  }
  return(data.frame(risk.score, risk.categ))
}

tcga.risk.score <- RiskEsti(exp.dat=tcga.lgg.exp, gene.set=active.k.vals.1se$entrez_id, risk.coef=active.k.vals.1se$coef, cut.off='Top')
tcga.lgg.cli.data <- merge(tcga.lgg.cli.data, tcga.risk.score, by='row.names')
colnames(tcga.lgg.cli.data)[1] <- 'patient_id'


# best threshold
hcPearsonWardAverCluster <- readRDS(file = '/result/section5/unFilter/hcPearsonWardAverCluster.rds')
lggCluster <- hcPearsonWardAverCluster[[5]]$consensusClass %>% as.data.frame() %>% rownames_to_column(var='patient_id')
colnames(lggCluster)[2] <- 'Clusters'
lggCluster <- lggCluster %>% mutate(Clusters = ifelse(Clusters == 4, 'Pos', 'Neg'))

tcga.lgg.cli.data <- merge(tcga.lgg.cli.data, lggCluster, by = 'patient_id')

DetermineOptThresPlot <- function(response, predictor, data, cutMethod = 'Youden'){
  
  rocRes <- roc(response=data[, response], predictor=data[, predictor])
  
  if(cutMethod == 'Youden'){
    # Maximization of Youden Index
    youden <- coords(rocRes, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
    
  }else{
    # Closest to the top left corner
    youden <- coords(rocRes, "best", best.method = "closest.topleft", ret = c("threshold", "sensitivity", "specificity"))
  }
  
  cat("Best threshold:", youden$threshold, '\n')
  
  # Draw the ROC curve and mark the optimal threshold
  ROCcurve <- ggroc(rocRes) + geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    geom_point(data = youden, aes(x = specificity, y = sensitivity), color = "red", size = 3) +
    annotate("text", x = youden$specificity - 0.15, y = youden$sensitivity - 0.05,
             label = paste("Threshold =", round(youden$threshold, 2))) +
    labs(x = "Specificity", y = "Sensitivity") + 
    theme_minimal() + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))
  
  return(ROCcurve)
}

# bestThresPlot <- DetermineOptThresPlot(response = 'Clusters', predictor = 'risk.score', data = tcga.lgg.cli.data, cutMethod = 'Youden')
bestThresPlot <- DetermineOptThresPlot(response = 'Clusters', predictor = 'risk.score', data = tcga.lgg.cli.data, cutMethod = 'closest.toplef')

tcga.lgg.cli.data$risk.categ <- ifelse(tcga.lgg.cli.data$risk.score > -3.29, 'C4-like', 'Others')


# > table(tcga.lgg.cli.data$Clusters, tcga.lgg.cli.data$risk.categ)
#     C4-like Others
# Neg      55    366
# Pos      73     12


riskScorePlot <- ggbarplot(data=tcga.lgg.cli.data, x='patient_id', y='risk.score', color=NA, fill='Clusters', sort.by.groups = FALSE, 
                xlab='Risk score for every patient', ylab='C4-like score', sort.val='desc') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7)) + 
  guides(fill=guide_legend(title=NULL)) + geom_vline(xintercept = 128, color = "red", linetype = "dashed")


ggsave(filename = '/result/section5/lgglike/cutoffChioce.pdf', 
       plot = bestThresPlot + riskScorePlot, width = 14, height = 7, units = 'cm')


saveRDS(tcga.lgg.cli.data, file='/data/LggRiskScores/tcga_lgg_risk_score.rds')


# survival plot 
source('/code/0.DataPreparation/survival_plot.R')
SurvivalPlot(survival.data=tcga.lgg.cli.data[, c('patient_id', 'os_time', 'os')], 
             sample.class=tcga.lgg.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lgg_os.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=tcga.lgg.cli.data[, c('patient_id', 'dss_time', 'dss')], 
             sample.class=tcga.lgg.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lgg_dss.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=tcga.lgg.cli.data[, c('patient_id', 'pfi_time', 'pfi')], 
             sample.class=tcga.lgg.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lgg_pfi.pdf', 
             out.file.path='/result/section5/lgglike/')

SurvivalPlot(survival.data=tcga.lgg.cli.data[, c('patient_id', 'dfi_time', 'dfi')], 
             sample.class=tcga.lgg.cli.data[, c('patient_id', 'risk.categ')], filename='tcga_lgg_dfi.pdf', 
             out.file.path='/result/section5/lgglike/')


################################uni_mul_cox
UnivariateCox <- function(cli.data, covariates){
  library('survival')
  #STEP1:构建单因素分析的对象
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time, status)~', x)));
  
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
Cox.function <- function(time, event, clinical.data, clinical.variate = NULL){
  ###设置工作环境
  options(stringsAsFactors = FALSE, warn = -1);
  suppressPackageStartupMessages(require(survival));
  
  ###判断协变量类型：数值型（num.covariate），非数值型(chara.covariate)。便于后续输出
  if(is.null(clinical.variate))
  {
    covariates  <- colnames(clinical.data)[-c(1:3)];
  }
  if(is.numeric(clinical.variate))
  {
    covariates  <- colnames(clinical.data)[clinical.variate];
  }
  num.variate <- NULL
  for(i in covariates)
  {
    if(is.numeric(clinical.data[, i]))
    {
      num.variate <- append(num.variate, i)
    }
  }
  chara.variate <- setdiff(covariates, num.variate)
  
  ####单因素cox分析函数：univariate.cox
  univariate.cox <- function(data, num, chara)
  {
    #STEP1:构建单因素分析的对象
    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(time, event)~', x)));
    
    #STEP2:单因素Cox分析
    univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)});
    
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
    );
    
    #STEP4:标准化输出格式
    for(i in num)
    {
      tmp <- univ_results[[i]];
      tmp$type <- " ";
      univ_results[[i]] <- tmp; 
    }
    
    for(i in chara)
    {
      tmp <- univ_results[[i]]
      tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
      tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
      tmp$variate <- tmp.variate;
      tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
      univ_results[[i]] <- tmp; 
    }
    
    #STEP5: list转化为数据框输出 	
    univ.result <- do.call(rbind.data.frame, univ_results);
    univ.result <- univ.result[,c(1,4,2,3)];
    colnames(univ.result) <- c('variate', 'type', 'univ HR (95% CI for HR)', 'univ p value')
    rownames(univ.result) <- NULL;
    return(univ.result)  
  };
  
  ####多因素cox分析函数：multivariate.cox
  multivariate.cox <- function(data, num, chara)
  {
    options(stringsAsFactors = FALSE);
    
    #STEP1:直接对所有选中的协变量进行多因素分析
    multiv_formula <- as.formula(paste("Surv(time, event)~", paste(covariates, collapse="+"), sep=""));
    multiv_model    <- coxph(multiv_formula, data = data);
    
    #STEP2:提取有用信息
    tmp <- summary(multiv_model);
    
    #提取p值，保留两位有效数字
    # p.value <- round(tmp$coefficients[ ,5], digits = 4);
    # p.value[which(p.value < 0.0001)] <- "<0.0001"; 
    p.value <- tmp$coefficients[ ,5]
    
    #提取beta值，这里得到的coefficients为矩阵，且有多行（每行对应一个协变量）
    #beta <- round(tmp$coefficients[ ,1], digits = 4);
    
    #提取风险比
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    
    #提取95%置信区间上下界
    HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4);
    HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"],digits = 4);
    
    #合并风险比HR和置信区间为一个内容
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
    
    variate <- rownames(tmp$coefficients);
    
    #整合输出内容为data.frame
    multiv_result <- as.data.frame(cbind(variate, HR, p.value));
    
    #STEP3:新建数据框储存多因素结果
    multiv.result <- NULL;
    
    for(i in num)
    {
      n.row <- grep(pattern = i, multiv_result$variate);
      tmp <- multiv_result[n.row, ];
      tmp$type <- " ";
      multiv.result <- rbind(multiv.result,tmp);
    }
    
    for(i in chara)
    {
      n.row <- grep(pattern = i, multiv_result$variate);
      tmp <- multiv_result[n.row, ];
      tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
      tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
      tmp$variate <- tmp.variate;
      tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
      multiv.result <- rbind(multiv.result,tmp);         
    }
    multiv.result <- multiv.result[,c(1,4,2,3)]
    colnames(multiv.result) <- c("variate", "type","multiv HR (95% CI for HR)", "multiv p value");
    rownames(multiv.result) <- NULL;
    
    return(multiv.result);
  };
  
  ###运行上述函数
  UniCoxPH <- univariate.cox(data = clinical.data, num = num.variate, chara = chara.variate);
  MultiCoxPH <- multivariate.cox(data = clinical.data, num = num.variate, chara = chara.variate);
  
  ###合并两个数据框
  cox.result <- merge(UniCoxPH, MultiCoxPH, by = c("variate", "type"), all = T);
  colnames(cox.result) <- c("variate", " ", "univ HR (95% CI for HR)", "univ p value","multiv HR (95% CI for HR)", "multiv p value");
  
  ###更改表格格式
  for(i in chara.variate)
  {
    tmp.row <- which(cox.result[,1] == i);
    tmp.vec <- c(i, rep(" ", times = 5));
    cox.result <- rbind(cox.result[1 : (tmp.row-1),], tmp.vec, cox.result[tmp.row : nrow(cox.result),]); 
  };
  cox.result[duplicated(cox.result[,1]),1] <- " ";
  rownames(cox.result) <- 1:nrow(cox.result);
  
  return(cox.result)
}


# CDKN2A/B Deletion, EGFR amplification
gene.cnv.alt <- read.table(file = '/data/GliomaData/all_thresholded.by_genes.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(gene.cnv.alt) <- gene.cnv.alt$Gene.Symbol
gene.cnv.alt <- gene.cnv.alt[, -c(1:3)]
colnames(gene.cnv.alt) <- gsub('\\.', '-', substr(colnames(gene.cnv.alt), 1, 15))
gene.cnv.alt <- as.data.frame(t(gene.cnv.alt))
gene.cnv.alt$CDKN2AB <- ifelse((gene.cnv.alt$CDKN2A == -2) | (gene.cnv.alt$CDKN2B == -2), 1, 0)
gene.cnv.alt$EGFR <- ifelse(gene.cnv.alt$EGFR == 2, 1, 0)


sub.cli.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
tcga_glioma_cli_mol <- readRDS(file = '/data/GliomaData/tcga_glioma_cli_mol.rds')
tcga_glioma_cli_mol <- tcga_glioma_cli_mol %>% mutate(bcr_patient_barcode = paste0(bcr_patient_barcode, '-01'))
tcgaPanCanSamples <- readRDS(file = '/data/tcgaPanCanSamples.rds')

sub.cli.data <- merge(sub.cli.data, tcgaPanCanSamples[, c('SAMPLE_BARCODE', 'SUBTYPE')], by.x = 'patient_id', by.y = 'SAMPLE_BARCODE')
sub.cli.data <- merge(sub.cli.data, tcga_glioma_cli_mol[, c('bcr_patient_barcode', 'MGMT_PROMOTER_STATUS')], 
                      by.x = 'patient_id', by.y = 'bcr_patient_barcode')

sub.cli.data$CDKN2AB <- gene.cnv.alt[sub.cli.data$patient_id, 'CDKN2AB']
sub.cli.data$EGFR <- gene.cnv.alt[sub.cli.data$patient_id, 'EGFR']

# > mean(sub.cli.data$age)
# [1] 42.99012
sub.cli.data <- sub.cli.data %>% mutate(age_categ = ifelse(age >= 43, '>= 43', '<43'), 
                                        gender = factor(gender,levels = c('FEMALE', 'MALE')), 
                                        # race = factor(race, levels = c('BLACK', 'OTHER', 'WHITE')),
                                        histological_type = factor(histological_type, levels = c('Oligodendroglioma', 'Oligoastrocytoma', 'Astrocytoma')), 
                                        histological_grade = factor(histological_grade, levels = c('G2', 'G3')),
                                        mgmt.status = factor(MGMT_PROMOTER_STATUS, levels = c('Unmethylated', 'Methylated')), 
                                        CDKN2AB = factor(CDKN2AB, levels = c(0, 1)), 
                                        EGFR = factor(EGFR, levels = c(0, 1)), 
                                        SUBTYPE = factor(SUBTYPE, levels = c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt')), 
                                        risk.categ = factor(risk.categ , levels = c('Others', 'C4-like')),
                                        Clusters = factor(Clusters, levels = c('Neg', 'Pos')))

cli.sig.char <- sub.cli.data %>% select(patient_id, os_time, os, age_categ, gender, histological_type, histological_grade, mgmt.status, CDKN2AB, EGFR, SUBTYPE, risk.categ) # race,
cli.sig.char <- cli.sig.char %>% rename(time = os_time, status = os)
# cli.sig.char <- sub.cli.data %>% select(patient_id, pfi_time, pfi, age_categ, gender, histological_type, histological_grade, mgmt.status, CDKN2AB, EGFR, SUBTYPE, risk.categ) # race,
# cli.sig.char <- cli.sig.char %>% rename(time = pfi_time, status = pfi)

cli.sig.char <- subset(cli.sig.char, !is.na(time) & !is.na(status))
uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[4:12])


# 多因素cox分析
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]
uni.mul.cox.res <- Cox.function(time=cli.sig.char$time, event=cli.sig.char$status, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:12))



