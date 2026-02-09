
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



# risk evaluation
load(file='/data/lgg_lasso_binomial_res.RData')
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

risk.score <- RiskEsti(exp.dat=exp.data, gene.set=active.k.vals.1se$entrez_id, risk.coef=active.k.vals.1se$coef, cut.off='Top')
sub.cli.data <- merge(sub.cli.data, risk.score, by='row.names')
colnames(sub.cli.data)[1] <- 'patient_id'


# survival plot 
source('/code/0.DataPreparation/survival_plot.R')
sub.cli.data <- sub.cli.data %>% mutate(pfs.event = ifelse(pfs.event == 'Yes', TRUE, FALSE))

SurvivalPlot(survival.data=sub.cli.data[, c('patient_id', 'pfs', 'pfs.event')], 
             sample.class=sub.cli.data[, c('patient_id', 'risk.categ')], filename='plgg_GSE107850_pfs.pdf', 
             out.file.path='/result/section5/lgglike/')


saveRDS(sub.cli.data, file='/data/LggRiskScores/plgg_GSE107850_risk_score.rds')



################################uni_mul_cox
UnivariateCox <- function(cli.data, covariates){
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
uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[4:11])


# 多因素cox分析
cli.sig.char <- subset(cli.sig.char, idh.status != 'undetermined')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$pfs, event=cli.sig.char$pfs.event, 
                                clinical.data=cli.sig.char, clinical.variate = c(4:11))




