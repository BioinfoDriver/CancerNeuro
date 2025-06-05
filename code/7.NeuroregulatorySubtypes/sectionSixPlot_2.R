
library('dplyr')

########
PlotSurv <- function(clinical.data, upper.time = NULL, xscale = 1, xlab = "Time", median.time = TRUE, 
                     surv.median.line = "none", HR = FALSE, risk.table = TRUE, pval = TRUE, 
                     conf.int = FALSE, main = NULL, ylab = "Survival probability") {
  
  #载入相关R包
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)
  
  #确定事件类型和时间的单位
  # survival.event <- survival.event[1];
  # unit.xlabel <- unit.xlabel[1];
  
  #如果设置upper.time，则去除生存时间超过upper.time的样本
  if (!is.null(upper.time)) clinical.data <- clinical.data[clinical.data$time <= upper.time,]
  
  # #选择日期格式 
  # xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
  # switch(unit.xlabel, year={xScale <- 365.25;}, month={xScale <- 30;}, week={xScale <- 7;}, day={xScale <- 1})
  # xLab <- xSL[which(xSL[,1]==xScale),2];
  
  #构造颜色
  if (!is.factor(clinical.data$sample.label)) 
    clinical.data$sample.label <- as.factor(clinical.data$sample.label)
  
  t.name <- levels(clinical.data$sample.label)
  
  if (length(t.name) > 6) stop("样本分组>6，超过函数接受范围")
  colors <- c('#C67326', '#896495', "#808080","#EA4335","#4285F4","#FBBC05","#34A853","#000000") # 顺序：灰，红，蓝，黄，绿，黑
  t.col <- colors[1:length(t.name)]
  
  # 构造生存对象
  km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
  
  # 计算HR值和95%CI
  if (length(t.name) == 2) {
    if (HR) {
      cox.obj <- coxph(Surv(time, event)~sample.label, data=clinical.data)
      tmp <- summary(cox.obj)
      HRs <- round(tmp$coefficients[ ,2], digits = 2)
      HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 2)
      HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 2)
      HRs <- paste0(HRs, " (", HR.confint.lower, "-", HR.confint.upper, ")")	  
    }
  }
  
  # 构造生存图像中图例显示文字
  legend.content <- substr(names(km.curves$strata), start = 14, stop = 1000)
  
  # x轴刻度单位转换
  if (is.numeric(xscale) | (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"))) {
    xscale = xscale
  } else {
    stop('xscale should be numeric or one of c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m").')
  }
  
  # 隐函数：转换生存时间单位
  .format_xticklabels <- function(labels, xscale){
    # 1 year = 365.25 days
    # 1 month = 365.25/12 = 30.4375 days
    if (is.numeric(xscale)) xtrans <- 1/xscale
    else
      xtrans <- switch(xscale,
                       d_m = 12/365.25,
                       d_y = 1/365.25,
                       m_d = 365.25/12,
                       m_y = 1/12,
                       y_d = 365.25,
                       y_m = 12,
                       1
      )
    round(labels*xtrans,2)
  }
  
  # 在图中添加中位生存时间及其95%CI,放在副标题位置
  subtitle <- NULL
  if (median.time) {
    if (is.numeric(xscale)) {
      median.km.obj = km.curves
    } else if (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m")) {
      clinical.data$time <- .format_xticklabels(labels = clinical.data$time, xscale = xscale)
      median.km.obj <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    }
    survival.time.info <- NULL
    survival.time.info <- rbind(survival.time.info, summary(median.km.obj)$table)
    median.survival <- round(survival.time.info[!duplicated(survival.time.info[,7:9]),7:9], digits = 2) # 注意：这里取得的置信区间上界可能为NA
    if (length(levels(clinical.data$sample.label)) == 1) {
      tmp1 <- levels(clinical.data$sample.label)
    } else {
      tmp1 <- do.call(rbind,strsplit(rownames(summary(median.km.obj)$table), split = "="))[,2]
    }
    tmp2 <- paste(median.survival[,1], "(", median.survival[,2], "-", median.survival[,3], ")")
    subtitle <- paste(tmp1, tmp2, sep = ":", collapse = "\n")
  }
  
  # ggsurvplot绘制生存图像
  ggsurv <- ggsurvplot(km.curves,               # survfit object with calculated statistics.
                       data = clinical.data,             # data used to fit survival curves.
                       palette = t.col,
                       
                       #图的主题构架
                       risk.table = risk.table,       # show risk table.
                       pval = pval,             # show p-value of log-rank test.
                       surv.median.line = surv.median.line,  # add the median survival pointer.
                       title = main,     #主标题名字
                       subtitle = subtitle, #副标题
                       font.main = 15,       #主标题字体大小              
                       xlab = xlab,   # customize X axis label.
                       ylab = ylab,   # customize Y axis label
                       xscale = xscale,
                       break.x.by = 365.25,
                       
                       
                       #图例设置
                       legend.title = "", #图例标题，一般不用，设为空
                       legend.labs = legend.content, #图例文字描述
                       legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
                       font.legend = 9,     #图例字体大小
                       
                       #risk table设置
                       tables.theme = theme_cleantable(),#table主题
                       risk.table.title = "No. at risk:",#table标题
                       risk.table.y.text.col = T, # 使用颜色代替Y轴文字
                       risk.table.y.text = FALSE, # Y轴不使用文字注释
                       tables.height = 0.15,      # table的高度
                       risk.table.fontsize = 3    # risk table内文字的大小
  )
  # HR标注的位置限定
  if (length(t.name) == 2) {
    if (HR) 
      ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = max(km.curves$time)/13, 
                                                     y = 0.15, size = 5, label = paste("HR=", HRs))
  } 
  #图的标题居中
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 10), 
                                     plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  #表的标题
  ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04),
                                       plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  
  # 判断分类的类数，如果只有两类，就不必计算两两之间的log rank p值
  if(length(t.name) > 2) {
    # 计算pairwise的log rank的p值
    res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data, p.adjust.method='none');
    pairwise.pvalue <- res$p.value;
    pairwise.pvalue[which(pairwise.pvalue < 0.001)] <- "<0.001";
    pairwise.pvalue[is.na(pairwise.pvalue)] <- "-"
    
    pairwise.pvalue <- apply(pairwise.pvalue, 2, substr, start = 1, stop = 6)
    # 添加表格
    tt <- ttheme_minimal(core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                         colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                         rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black"))
    )
    pairwise.table <- tableGrob(pairwise.pvalue, theme = tt)
    ggsurv <- ggarrange(ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5)),
                        pairwise.table, nrow=2, heights = c(2,0.5),
                        labels = c("","pairwise comparisons"),
                        hjust = 0, font.label = list(size = 15, face = "plain"))
  }
  
  ggsurv	
}

########


setwd('/result/section6/consensusCluster_New/New3')

panCanCluster <- read.csv(file='panCanClusterNew3.k=5.consensusClass.csv', header = F, sep = ',', stringsAsFactors = FALSE)
panCanCluster <- panCanCluster %>% rename(SAMPLE_BARCODE = V1, Clusters = V2)  %>% mutate(patientID = substr(SAMPLE_BARCODE, 1, 12))

panCanCliData <- readRDS(file = '/data/tcgaPanCanCliData.rds')
panCanCliData <- panCanCliData %>% tibble::rownames_to_column(var = "patientID")

panCanCliData <- panCanCliData %>% right_join(panCanCluster, by = 'patientID')


osCanTypes <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
                'LIHC','LUAD','LUSC','MESO','OV','PAAD','PRAD','READ','SARC','SKCM','STAD','THCA','UCEC','UCS','UVM')
panCanCliData <- subset(panCanCliData, cancer_type %in% osCanTypes & !is.na(os) & !is.na(os_time))

###
cliData <- panCanCliData[, c('patientID', 'os', 'os_time', 'Clusters')]
colnames(cliData) <- c('patientID', 'event', 'time', 'sample.label')


pdf(file = '/result/section6/consensusCluster_New/New3/panCan5clusterSurplot.pdf')

surplot <- try(PlotSurv(clinical.data = cliData, xscale = 'd_y', xlab = 'Time in years', 
                        ylab = 'Overall survival', risk.table = FALSE, median.time = FALSE))

print(surplot)
dev.off()

###
panCanCliData <- split.data.frame(panCanCliData, f = panCanCliData$cancer_type)
clusterSurplot <- lapply(panCanCliData, function(cliData){
  
  disease <- cliData$cancer_type[1]
  
  cliData <- cliData[, c('patientID', 'os', 'os_time', 'Clusters')]
  colnames(cliData) <- c('patientID', 'event', 'time', 'sample.label')
  
  
  surplot <- try(PlotSurv(clinical.data = cliData, xscale = 'd_y', xlab = 'Time in years', 
                          ylab = 'Overall survival', risk.table = FALSE, median.time = FALSE, main = disease))
  
  return(surplot)
})



pdf(file = '/result/section6/consensusCluster_New/New3/5clusterSurplot.pdf')


for(i in seq(length(clusterSurplot))){
  print(clusterSurplot[[i]])
}

dev.off()

