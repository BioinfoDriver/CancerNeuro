
library('DESeq2')
library('dplyr')

# read count data
exp.count <- read.csv(file = 'D:/CancerNeuroscience/data/GliomaData/CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt', 
                      header = T, sep = '\t', stringsAsFactors = F)

exp.count <- exp.count %>% tibble::column_to_rownames(var = 'gene_name')

# read clinical data
cli.data <- read.csv(file = 'D:/CancerNeuroscience/data/GliomaData/CGGA.mRNAseq_325_clinical.20200506.txt', 
                     header = T, sep = '\t', stringsAsFactors = F)

cli.data <- cli.data %>% tibble::column_to_rownames(var = 'CGGA_ID')
cli.data <- cli.data[colnames(exp.count), ]

cli.data <- cli.data %>% mutate(disType = case_when(PRS_type == 'Primary' & Histology == 'GBM' ~ 'pGBM', 
                                                    PRS_type == 'Recurrent' & Histology == 'rGBM' ~ 'rGBM',
                                                    PRS_type == 'Primary' & Histology %in% c('A', 'AA', 'AO', 'AOA', 'O', 'OA') ~ 'pLGG',
                                                    PRS_type == 'Recurrent' & Histology %in% c('rA', 'rAA', 'rAO', 'rAOA', 'rO', 'rOA') ~ 'rLGG',
                                                    PRS_type == 'Secondary' & Histology == 'sGBM' ~ 'sGBM'))

# The variance stabilizing transformation
sam.info <- cli.data[, 'disType', FALSE]

dds <- DESeqDataSetFromMatrix(countData = exp.count, colData=sam.info, design = ~ 1)
vsd <- vst(dds, blind = FALSE)
norm.exp.data <- assay(vsd)


sub.cli.data <- subset(cli.data, disType == 'pLGG')
exp.data <- norm.exp.data[, intersect(colnames(norm.exp.data), rownames(sub.cli.data))]
sub.cli.data <- sub.cli.data[intersect(colnames(exp.data), rownames(sub.cli.data)), ]


# risk evaluation
load(file='D:/CancerNeuroscience/Github/data/lgg_lasso_binomial_res.RData')
RiskEsti <- function(exp.dat, gene.set, risk.coef, cut.off=NULL){
  
  # exp.dat <- as.matrix(exp.dat[intersect(rownames(exp.dat), gene.set), ])
  exp.dat <- as.matrix(exp.dat[gene.set, ])
  
  risk.score <- crossprod(exp.dat, matrix(risk.coef, nrow=length(risk.coef)))[, 1]
  if(!is.null(cut.off)){
    risk.categ <- ifelse(risk.score >= quantile(risk.score, 0.75), 'high risk', 'low risk')
    
  }else{
    risk.categ <- ifelse(risk.score >= median(risk.score), 'high risk', 'low risk')
    
  }
  return(data.frame(risk.score, risk.categ))
}

risk.score <- RiskEsti(exp.dat=exp.data, gene.set=active.k.vals.1se$symbol, risk.coef=active.k.vals.1se$coef, cut.off='Top')
sub.cli.data <- merge(sub.cli.data, risk.score, by='row.names')
colnames(sub.cli.data)[1] <- 'patient_id'


# survival plot 
source('D:/CancerNeuroscience/Github/code/0.DataPreparation/survival_plot.R')

SurvivalPlot(survival.data=sub.cli.data[, c('patient_id', 'OS', 'Censor..alive.0..dead.1.')], 
             sample.class=sub.cli.data[, c('patient_id', 'risk.categ')], filename='plgg_325_os.pdf', 
             out.file.path='D:/CancerNeuroscience/Github/result/section5/lgglike/')

saveRDS(sub.cli.data, file='D:/CancerNeuroscience/Github/data/LggRiskScores/plgg_325_risk_score.rds')



################################uni_mul_cox
UnivariateCox <- function(cli.data, covariates){
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

sub.cli.data <- readRDS(file='D:/CancerNeuroscience/Github/data/LggRiskScores/plgg_325_risk_score.rds')

# > mean(sub.cli.data$Age, na.rm = T)
# [1] 40.66667
sub.cli.data <- sub.cli.data %>% mutate(os = Censor..alive.0..dead.1., os_time = OS, 
                                        age_categ = ifelse(Age >= 41, '>= 41', '<41'),
                                        gender = factor(Gender,levels = c('Female', 'Male')), 
                                        Histology = factor(Histology, levels = c('O', 'OA', 'AO', 'AOA', 'A', 'AA')), 
                                        Grade = factor(Grade, levels = c('WHO II', 'WHO III')),
                                        idh.status = factor(IDH_mutation_status, levels = c('Mutant', 'Wildtype')), 
                                        X1p19q.codeletion = factor(X1p19q_codeletion_status, levels = c('Codel', 'Non-codel')),
                                        mgmt.status = factor(MGMTp_methylation_status, levels = c('un-methylated', 'methylated')), 
                                        TMZ.treated = factor(Chemo_status..TMZ.treated.1.un.treated.0., levels = c(0, 1)),
                                        Radio.status = factor(Radio_status..treated.1.un.treated.0., levels = c(0, 1)),
                                        risk.categ = factor(risk.categ , levels = c('low risk', 'high risk')))

cli.sig.char <- sub.cli.data %>% select(patient_id, os_time, os, age_categ, gender, Grade, idh.status, X1p19q.codeletion, 
                                        risk.categ, mgmt.status, Histology) # , TMZ.treated, Radio.status

cli.sig.char <- subset(cli.sig.char, !is.na(os_time) & !is.na(os))
uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[4:11])


# 多因素cox分析
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:11))

