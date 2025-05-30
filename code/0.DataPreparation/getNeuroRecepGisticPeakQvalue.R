
library('magrittr')
library('GenomicRanges')

setwd('/data/Firehose')

###########################
# fileList <- dir()
# renameFile <- gsub('-', '_', gsub('.CopyNumber', '_Gistic2.tar.gz', do.call(rbind, strsplit(fileList, '_'))[, 2]))
# 
# 
# for(i in seq(length(fileList))){
#   file.rename(fileList[i], renameFile[i])
# }
# 
# ###########################
# fileList <- dir()
# 
# for(i in seq(length(fileList))){
#   
#   untar(tarfile = fileList[i], exdir = gsub('.tar.gz', '', fileList[i]))
# 
# }
###########################

fileNames <- dir(recursive = TRUE, pattern = 'all_lesions.conf_99.txt$')

gistic2Peak <- lapply(fileNames, function(fileName){
  
  peakRegion <- read.csv(file = fileName, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  peakRegion <- peakRegion[1:(nrow(peakRegion)/2), 1:9]
  
  peakRegion <- peakRegion %>% mutate(Descriptor = trimws(Descriptor), Wide.Peak.Limits = trimws(Wide.Peak.Limits))
  peakRegion$Wide.Peak.Limits <- gsub("\\(probes\\s\\d+:\\d+\\)", "", peakRegion$Wide.Peak.Limits)
  peakRegion$Disease <- sub("_.*_.*", "", fileName)

  peakRegion <- peakRegion %>% mutate(chr = sub(":.*", "", Wide.Peak.Limits), 
                                      start = as.numeric(sub(".*:", "", sub("-.*", "", Wide.Peak.Limits))), 
                                      end = as.numeric(sub(".*-", "", Wide.Peak.Limits)))
  return(peakRegion)
})

gistic2Peak <- do.call(rbind, gistic2Peak)
gistic2Peak <- subset(gistic2Peak, q.values < 0.25)

gistic2Peak <- makeGRangesFromDataFrame(gistic2Peak, keep.extra.columns=TRUE)


###########################
nrReceptors <- readRDS(file = '/data/neurotransmitterReceptors.rds')
geneInfo <- readRDS(file = '/data/geneInfo.rds')


otherGenes <- subset(geneInfo, !(gene_name %in% nrReceptors$Approved.symbol))
nrReceptors <- merge(nrReceptors, select(geneInfo, gene_name, seqname, start, end), 
                     by.x = 'Approved.symbol', by.y = 'gene_name')

nrReceptors <- makeGRangesFromDataFrame(nrReceptors, seqnames.field = 'seqname', keep.extra.columns=TRUE)
otherGenes <- makeGRangesFromDataFrame(otherGenes, seqnames.field = 'seqname', keep.extra.columns=TRUE)

###########################
rgOverlaps  <- findOverlaps(nrReceptors, gistic2Peak, type = 'within', select = 'all')


ovRegions <- nrReceptors[queryHits(rgOverlaps)]
attributes <- mcols(gistic2Peak)[subjectHits(rgOverlaps), ]

mcols(ovRegions) <- cbind.data.frame(mcols(ovRegions), attributes)


nrGisticPeakQvalue <- mcols(ovRegions)
nrGisticPeakQvalue$Direction <- sub("\\s.*", "", nrGisticPeakQvalue$Unique.Name)

saveRDS(nrGisticPeakQvalue, file = '/data/nrGisticPeakQvalue.rds')


###########################
ogOverlaps  <- findOverlaps(otherGenes, gistic2Peak, type = 'within', select = 'all')


ovRegions <- otherGenes[queryHits(ogOverlaps)]
attributes <- mcols(gistic2Peak)[subjectHits(ogOverlaps), ]

mcols(ovRegions) <- cbind.data.frame(mcols(ovRegions), attributes)


ogGisticPeakQvalue <- mcols(ovRegions)
ogGisticPeakQvalue$Direction <- sub("\\s.*", "", ogGisticPeakQvalue$Unique.Name)

saveRDS(ogGisticPeakQvalue, file = '/data/ogGisticPeakQvalue.rds')






