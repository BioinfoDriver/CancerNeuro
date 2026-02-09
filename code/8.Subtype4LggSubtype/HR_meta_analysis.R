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

uniCox <- data.frame(dataset=c('TCGA RNAseq', 'CGGA RNAseq-1', 'CGGA RNAseq-2', 'CGGA array', 'GSE55918'), 
                        HR=c(4.86, 3.06, 4.75, 4.34, 2.47), 
                        lower=c(3.36, 2.03, 2.85, 2.63, 1.86), upper=c(7.03, 4.60, 7.91, 7.18, 3.29), 
                        pvalue=c('4.55e-17', '8.29e-08', '2.19e-09', '1.01e-08', '5.45e-10'))



# RmetaHrMetaAnalysis(uniCox, 'fixed', TRUE, 'D:/CancerNeuroscience/Github/result/section5/lgglike/metaUniHRmetafor.pdf')
MetaforHrMetaAnalysis(uniCox, 'FE', TRUE, 'D:/CancerNeuroscience/Github/result/section5/lgglike/metaUniHRrmeta.pdf')

# hr        ci.lb        ci.ub         pval 
# 3.432082e+00 2.886915e+00 4.080200e+00 2.291211e-44


mulCox <- data.frame(dataset=c('TCGA RNAseq', 'CGGA RNAseq-1', 'CGGA RNAseq-2', 'CGGA array', 'GSE55918'), 
                      HR=c(1.81, 2.55, 2.02, 4.62, 2.84), 
                      lower=c(1.03, 1.43, 1.07, 2.37, 2.10), upper=c(3.18, 4.56, 3.84, 9.04, 3.84), 
                      pvalue=c('3.83e-02', '1.52e-03', '3.05e-02', '7.57e-06', '1.63e-11'))

# RmetaHrMetaAnalysis(mulCox, 'fixed', TRUE, 'D:/CancerNeuroscience/Github/result/section5/lgglike/metaMulHRmetafor.pdf')
MetaforHrMetaAnalysis(mulCox, 'FE', TRUE, 'D:/CancerNeuroscience/Github/result/section5/lgglike/metaMulHRrmeta.pdf')

# hr        ci.lb        ci.ub         pval 
# 2.652620e+00 2.141059e+00 3.286408e+00 4.479135e-19 