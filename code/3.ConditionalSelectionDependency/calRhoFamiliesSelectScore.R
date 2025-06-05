
library(select)
#######################
load(file = '/data/dataForSelectCal.RData')

# Run select
alpi = select::select(M=gam, sample.class = sample.class, alteration.class = alteration.class, 
                      folder='/result/section2/', 
                      min.feature.support = 8, min.feature.freq = 0.001, 
                      n.permut = 1000, r.seed = 110, n.cores=32, calculate_APC_threshold = FALSE)

#######################

TP.FP.stats <- function (A, res = 1000) 
{
  A = A[which(A$APC > 1e-08), , drop = FALSE]
  if (nrow(A) < 1) 
    return(NULL)
  if (nrow(A[A$FDR, ]) < 1) 
    return(NULL)
  if (nrow(A[!A$FDR, ]) < 1) 
    return(NULL)
  t1 = min(A$APC) + ((0:res) * ((max(A$APC) - min(A$APC))/res))
  t2 = unique(A[A$FDR, , drop = FALSE]$APC)
  t = sort(unique(c(t1, t2)))
  AA = data.frame(t = t)
  AA$TN = ecdf(A[!A$FDR, , drop = FALSE]$APC)(AA$t)
  AA$FN = ecdf(A[A$FDR, , drop = FALSE]$APC)(AA$t)
  AA$TP = 1 - AA$FN
  AA$FP = 1 - AA$TN
  AA$F1 = 2 * AA$TP/(2 * AA$TP + AA$FN + AA$FP)
  AA$MCC = ((AA$TP * AA$TN) - (AA$FP * AA$FN))/sqrt((AA$TP + 
                                                       AA$FP) * (AA$TP + AA$FN) * (AA$TN + AA$FP) * (AA$TN + 
                                                                                                       AA$FN))
  return(as.vector(AA))
}


TP.FP.stat <- function (A, th) 
{
  A = A[which(A$APC > 1e-08), , drop = FALSE]
  t = th
  TN = ecdf(A[!A$FDR, , drop = FALSE]$APC)(t)
  FN = ecdf(A[A$FDR, , drop = FALSE]$APC)(t)
  TP = 1 - FN
  FP = 1 - TN
  F1 = 2 * TP/(2 * TP + FN + FP)
  MCC = ((TP * TN) - (FP * FN))/sqrt((TP + FP) * (TP + FN) * 
                                       (TN + FP) * (TN + FN))
  return(c(t = t, TN = TN, FN = FN, TP = TP, FP = FP, F1 = F1, 
           MCC = MCC))
}


get_thresh.2 <- function (x, A) 
{
  TP.FP.thresh_v1 = x[x[, "FP"] <= 0.05, , drop = FALSE]
  TP.FP.thresh_v2 = x[x[, "FP"] >= 0.05, , drop = FALSE]
  if (is.null(x)) 
    return(NA)
  if (is.null(dim(x))) 
    return(NA)
  if (nrow(x) == 0) 
    return(NA)
  if (nrow(TP.FP.thresh_v1) == 0) 
    return(x[nrow(x), ])
  if (nrow(TP.FP.thresh_v2) == 0) 
    return(TP.FP.thresh_v1[1, ])
  th = (TP.FP.thresh_v1[1, "t"] + TP.FP.thresh_v2[nrow(TP.FP.thresh_v2), 
                                                  "t"])/2
  th = (TP.FP.thresh_v1[1, "t"])
  th = TP.FP.thresh_v2[nrow(TP.FP.thresh_v2), "t"]
  return(TP.FP.stat(A, th))
}


establish_APC_threshold <- function(alpi) 
{
  if (nrow(alpi) < 25) {
    stats = NULL
    thresh = rep(NA, 7)
  }
  else {
    stats = TP.FP.stats(alpi, 10000)
    stats = as.data.frame(stats)
    
    thresh = get_thresh.2(stats, alpi)
  }
  return(list(stats = stats, thresh = thresh))
}


#########
TP.FP.thr <- establish_APC_threshold(alpi)
threshold <- TP.FP.thr$thresh["t"]
alpi$APC_good <- FALSE
alpi$APC_good[which(alpi$APC >= threshold)] <- TRUE

#########
saveRDS(alpi, file = '/data/nRSelectScoreRes.rds')


