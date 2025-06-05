# Function to fit meta-analytic equal-, fixed-, and random-effects models 
# https://stats.stackexchange.com/questions/343316/hazard-ratio-meta-analysis
# https://www.bmj.com/content/343/bmj.d2304
MetaforHrMetaAnalysis <- function(dat, method=c('EE', 'FE', 'REML'), ifplot, fname){
  method <- match.arg(method)
  
  library('metafor')
  library('dplyr')
  
  colnames(dat) <- c('study', 'hr', 'ci.lb', 'ci.ub', 'pval')
  dat <- dat %>% mutate(yi = log(hr), sei = (log(ci.ub)-log(ci.lb))/(2*1.96))
  
  # whether an equal-, a fixed or random-effects model should be fitted.
  res <- rma(yi=yi, sei=sei, data=dat, method=method, slab=study)
  
  meta.sum <- exp(c(res$b[1, 1], res$ci.lb, res$ci.ub))
  meta.sum <- c(meta.sum, res$pval)
  names(meta.sum) <- c('hr', 'ci.lb', 'ci.ub', 'pval')
  
  # plot
  if(ifplot){
    pdf(fname)
    forest(x=res, annotate=TRUE, header=c('Author(s) and Year', 'HR [95% CI]'), refline=1, 
           xlab='Hazard ratio', mlab='Overall', ilab=dat$pval, ilab.xpos=7, ilab.pos=2,
           colout='#3182bd', col='#a50f15', transf=exp)
    dev.off()
  }
  
  return(meta.sum)
}


RmetaHrMetaAnalysis <- function(dat, method=c("fixed", "random"), ifplot, fname){
  method <- match.arg(method)
  
  library('rmeta')
  library('dplyr')
  
  colnames(dat) <- c('study', 'hr', 'ci.lb', 'ci.ub', 'pval')
  dat <- dat %>% mutate(yi = log(hr), sei = (log(ci.ub)-log(ci.lb))/(2*1.96))
  
  # whether a fixed- or random-effects model should be fitted.
  res <- meta.summaries(d = yi, se = sei, method = method, logscale=TRUE, names=study, data=dat)
  
  meta.sum <- exp(c(res$summary, res$summary - 1.96 * res$se.summary, res$summary + 1.96 * res$se.summary))
  meta.sum <- c(meta.sum, res$test[2])
  names(meta.sum) <- c('hr', 'ci.lb', 'ci.ub', 'pval')
  
  # plot
  if(ifplot){
    pdf(fname)
    metaplot(mn=dat$yi, se=dat$sei, labels=dat$study, xlab='Hazard Ratio', 
             summn = res$summary, sumse = res$se.summary, sumnn= 1/res$se.summary^2, xlim=c(-1, 3), summlabel="Overall",
             zero=0, colors=meta.colors(box="#3182bd",lines="#a50f15", zero="red", summary="black",text="black"), xaxt='n')
    axis(1, at=log(c(0.5,1,2,4,8,16)), labels=c(0.5,1,2,4,8,16))
    dev.off()
  }
  
  return(meta.sum)
}

exam.data <- data.frame(dataset=c('TCGA RNAseq', 'CGGA RNAseq-1', 'CGGA RNAseq-2', 'CGGA array', 'GSE55918'), 
                        HR=c(6.46, 3.44, 4.96, 4.24, 2.67), 
                        lower=c(4.38, 2.27, 2.90, 2.53, 1.95), upper=c(9.53, 5.22, 8.47, 7.12, 3.66), 
                        pvalue=c('4.61e-21', '6.45e-09', '4.62e-09', '4.64e-08', '1.15e-09'))

RmetaHrMetaAnalysis(exam.data, 'fixed', TRUE, '/result/section6/lgglike/metaHRmetafor.pdf')
MetaforHrMetaAnalysis(exam.data, 'FE', TRUE, '/result/section6/lgglike/metaHRrmeta.pdf')

# hr        ci.lb        ci.ub         pval 
# 3.892472e+00 3.239890e+00 4.676498e+00 9.658006e-48 