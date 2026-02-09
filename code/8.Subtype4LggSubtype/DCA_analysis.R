
install.packages("dcurves")
library('dcurves')
library('survival')


# TCGA-LGG
tcga.lgg.cli.data <- readRDS(file='/data/LggRiskScores/tcga_lgg_risk_score.rds')
tcga.lgg.cli.data <- tcga.lgg.cli.data %>% select(os_time, os, risk.categ, age, histological_grade) %>% na.omit()

# Run the cox model
coxC4 <- coxph(Surv(os_time, os) ~ risk.categ, data = tcga.lgg.cli.data)
coxAge <- coxph(Surv(os_time, os) ~ age, data = tcga.lgg.cli.data)
coxGrade <- coxph(Surv(os_time, os) ~ histological_grade, data = tcga.lgg.cli.data)
coxFull <- coxph(Surv(os_time, os) ~ risk.categ + age + histological_grade, data = tcga.lgg.cli.data)


tcga.lgg.cli.data <-tcga.lgg.cli.data %>% 
  mutate(os_C4 = 1 - summary(survfit(coxC4, newdata = tcga.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Age = 1 - summary(survfit(coxAge, newdata = tcga.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Grade = 1 - summary(survfit(coxGrade, newdata = tcga.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Full = 1 - summary(survfit(coxFull, newdata = tcga.lgg.cli.data), times = 365*3)$surv[1, ])

tcgaDCA <- dca(Surv(os_time, os) ~ os_C4 + os_Age + os_Grade + os_Full, data = tcga.lgg.cli.data,
    time = 365*3, thresholds = seq(0, 0.99, by = 0.01),
    label = list(os_C4 = "C4-like", os_Age = "Age", os_Grade = "Grade", os_Full = 'C4-like/Age/Grade')) %>% 
  plot(smooth = TRUE)  + labs(title = 'TCGA-LGG cohort') + 
  theme(legend.position = c(0.9, 0.8), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                              legend.title = element_text(size = 7), legend.text = element_text(size = 6))


# CGGA-mRNAseq1 cohort
mseq1.lgg.cli.data <- readRDS(file='/data/LggRiskScores/plgg_693_risk_score.rds')
mseq1.lgg.cli.data <- mseq1.lgg.cli.data %>% select(Grade, Age, OS, Censor..alive.0..dead.1., risk.categ) %>% 
  rename(os_time = OS, os = Censor..alive.0..dead.1.) %>% na.omit()

# Run the cox model
coxC4 <- coxph(Surv(os_time, os) ~ risk.categ, data = mseq1.lgg.cli.data)
coxAge <- coxph(Surv(os_time, os) ~ Age, data = mseq1.lgg.cli.data)
coxGrade <- coxph(Surv(os_time, os) ~ Grade, data = mseq1.lgg.cli.data)
coxFull <- coxph(Surv(os_time, os) ~ risk.categ + Age + Grade, data = mseq1.lgg.cli.data)


mseq1.lgg.cli.data <-mseq1.lgg.cli.data %>% 
  mutate(os_C4 = 1 - summary(survfit(coxC4, newdata = mseq1.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Age = 1 - summary(survfit(coxAge, newdata = mseq1.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Grade = 1 - summary(survfit(coxGrade, newdata = mseq1.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Full = 1 - summary(survfit(coxFull, newdata = mseq1.lgg.cli.data), times = 365*3)$surv[1, ])

mRNAseq1DCA <- dca(Surv(os_time, os) ~ os_C4 + os_Age + os_Grade + os_Full, data = mseq1.lgg.cli.data,
               time = 365*3, thresholds = seq(0, 0.60, by = 0.01),
               label = list(os_C4 = "C4-like", os_Age = "Age", os_Grade = "Grade", os_Full = 'C4-like/Age/Grade')) %>% 
  plot(smooth = TRUE) + labs(title = 'CGGA-mRNAseq1 cohort') + 
  theme(legend.position = c(0.9, 0.8), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
                              legend.title = element_text(size = 7), legend.text = element_text(size = 6))


# CGGA-mRNAseq2 cohort
mseq2.lgg.cli.data <- readRDS(file='/data/LggRiskScores/plgg_325_risk_score.rds')
mseq2.lgg.cli.data <- mseq2.lgg.cli.data %>% select(Grade, Age, OS, Censor..alive.0..dead.1., risk.categ) %>% 
  rename(os_time = OS, os = Censor..alive.0..dead.1.) %>% na.omit()

# Run the cox model
coxC4 <- coxph(Surv(os_time, os) ~ risk.categ, data = mseq2.lgg.cli.data)
coxAge <- coxph(Surv(os_time, os) ~ Age, data = mseq2.lgg.cli.data)
coxGrade <- coxph(Surv(os_time, os) ~ Grade, data = mseq2.lgg.cli.data)
coxFull <- coxph(Surv(os_time, os) ~ risk.categ + Age + Grade, data = mseq2.lgg.cli.data)


mseq2.lgg.cli.data <-mseq2.lgg.cli.data %>% 
  mutate(os_C4 = 1 - summary(survfit(coxC4, newdata = mseq2.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Age = 1 - summary(survfit(coxAge, newdata = mseq2.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Grade = 1 - summary(survfit(coxGrade, newdata = mseq2.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Full = 1 - summary(survfit(coxFull, newdata = mseq2.lgg.cli.data), times = 365*3)$surv[1, ])

mRNAseq2DCA <- dca(Surv(os_time, os) ~ os_C4 + os_Age + os_Grade + os_Full, data = mseq2.lgg.cli.data,
                   time = 365*3, thresholds = seq(0, 0.99, by = 0.01),
                   label = list(os_C4 = "C4-like", os_Age = "Age", os_Grade = "Grade", os_Full = 'C4-like/Age/Grade')) %>% 
  plot(smooth = TRUE) + labs(title = 'CGGA-mRNAseq2 cohort') + 
  theme(legend.position = c(0.9, 0.8), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))


# CGGA-mRNA-array1 cohort
arr1.lgg.cli.data <- readRDS(file='/data/LggRiskScores/plgg_301_risk_score.rds')
arr1.lgg.cli.data <- arr1.lgg.cli.data %>% select(Grade, Age, OS, Censor..alive.0..dead.1., risk.categ) %>% 
  rename(os_time = OS, os = Censor..alive.0..dead.1.) %>% na.omit()

# Run the cox model
coxC4 <- coxph(Surv(os_time, os) ~ risk.categ, data = arr1.lgg.cli.data)
coxAge <- coxph(Surv(os_time, os) ~ Age, data = arr1.lgg.cli.data)
coxGrade <- coxph(Surv(os_time, os) ~ Grade, data = arr1.lgg.cli.data)
coxFull <- coxph(Surv(os_time, os) ~ risk.categ + Age + Grade, data = arr1.lgg.cli.data)


arr1.lgg.cli.data <-arr1.lgg.cli.data %>% 
  mutate(os_C4 = 1 - summary(survfit(coxC4, newdata = arr1.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Age = 1 - summary(survfit(coxAge, newdata = arr1.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Grade = 1 - summary(survfit(coxGrade, newdata = arr1.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Full = 1 - summary(survfit(coxFull, newdata = arr1.lgg.cli.data), times = 365*3)$surv[1, ])

array1DCA <- dca(Surv(os_time, os) ~ os_C4 + os_Age + os_Grade + os_Full, data = arr1.lgg.cli.data,
                 time = 365*3, thresholds = seq(0, 0.99, by = 0.01),
                 label = list(os_C4 = "C4-like", os_Age = "Age", os_Grade = "Grade", os_Full = 'C4-like/Age/Grade')) %>% 
  plot(smooth = TRUE) + labs(title = 'CGGA-mRNA-array1 cohort') + 
  theme(legend.position = c(0.9, 0.8), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))


# mRNA-array2 cohort
arr2.lgg.cli.data <- readRDS(file='/data/LggRiskScores/plgg_GSE55918_risk_score.rds')
arr2.lgg.cli.data <- arr2.lgg.cli.data %>% select(Tumor.grade, Age, Survival.time.months, Censored, risk.categ) %>% 
  rename(Grade = Tumor.grade, os_time = Survival.time.months, os = Censored) %>% na.omit() %>% 
  mutate(Grade = factor(Grade, levels = c('G2', 'G3')))

# Run the cox model
coxC4 <- coxph(Surv(os_time, os) ~ risk.categ, data = arr2.lgg.cli.data)
coxAge <- coxph(Surv(os_time, os) ~ Age, data = arr2.lgg.cli.data)
coxGrade <- coxph(Surv(os_time, os) ~ Grade, data = arr2.lgg.cli.data)
coxFull <- coxph(Surv(os_time, os) ~ risk.categ + Age + Grade, data = arr2.lgg.cli.data)


arr2.lgg.cli.data <-arr2.lgg.cli.data %>% 
  mutate(os_C4 = 1 - summary(survfit(coxC4, newdata = arr2.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Age = 1 - summary(survfit(coxAge, newdata = arr2.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Grade = 1 - summary(survfit(coxGrade, newdata = arr2.lgg.cli.data), times = 365*3)$surv[1, ],
         os_Full = 1 - summary(survfit(coxFull, newdata = arr2.lgg.cli.data), times = 365*3)$surv[1, ])

array2DCA <- dca(Surv(os_time, os) ~ os_C4 + os_Age + os_Grade + os_Full, data = arr2.lgg.cli.data,
                 time = 365*3, thresholds = seq(0, 0.99, by = 0.01),
                 label = list(os_C4 = "C4-like", os_Age = "Age", os_Grade = "Grade", os_Full = 'C4-like/Age/Grade')) %>% 
  plot(smooth = TRUE) + labs(title = 'GEO-mRNA-array2 cohort') + 
  theme(legend.position = c(0.9, 0.8), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))


ggsave(plot = patchwork::wrap_plots(tcgaDCA,mRNAseq1DCA,mRNAseq2DCA,array1DCA,array2DCA, ncol = 3), 
       width = 20, height = 12, units = 'cm',
       filename = '/result/section5/lgglike/DCA_analysis.pdf')


