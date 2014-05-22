
library(stringr)
library(DESeq)
library(plyr)
library(reshape2)
library(ggplot2)
Stangle('s_functions.Rnw');source('s_functions.R')

load('./intermediate_results/linccounts.RData')
load('./intermediate_results/fullTable.RData')

lincsunion <- subset(linccounts, transcriptvar == 'merged')
row.names(lincsunion) <- lincsunion$GeneSymbol


# I would like to add these data to the Full data set and reanalyze with them in. 
# That will make for a fairer comparison.

datacolumns <- c('A1','A2','A3','A4','A5',
                 'C1','C2','C3','C4','C5')
conds <- c('ALL','GFP','PLUSS','NEG','DIFF',
           'ALL','GFP','PLUSS','NEG','DIFF')

lincsfordeseq <- rbind( lincsunion[, datacolumns], fullTable[, datacolumns] ) 



cdslincs <- newCountDataSet( lincsfordeseq, conds) 
cdslincs <- estimateSizeFactors(cdslincs)
cdslincs <- estimateDispersions( cdslincs , method='per-condition',
	sharingMode='maximum', fitType='parametric')



f.deseq_do_custom <- function(countdataset, condpair, 
	basename = 'nbinomresult') {
    # The function could create trouble if condition names have periods.
    message(paste(basename,':','Calculating',condpair[1],condpair[2]))
    nbinomresult<-nbinomTest(countdataset, condpair[1], condpair[2])
    filename<-paste( paste0('./intermediate_results/',basename),
                             condpair[1] , condpair[2], sep = '.')
    write.table(nbinomresult, file=filename, quote=FALSE, sep='\t')
    return( list(condA=condpair[1], condB=condpair[2],
                 filename=filename))
}



nbinomindex.wlincs <- apply(X = combn(levels(conditions(cdslincs)),m=2),
                   MARGIN=2, FUN= function (x) f.deseq_do_custom(cdslincs,x,'nbinomwlincs'))
save(file='./intermediate_results/nbinomindex.wlincs.RData', list='nbinomindex.wlincs')





