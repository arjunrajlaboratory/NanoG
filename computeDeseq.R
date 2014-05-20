require(DESeq)

f.deseq_do_andwrite <- function(countdataset, condpair,
                                fileprefix ) {
    message( paste( 'Testing', condpair[1], 'vs', condpair[2]))
    nbinomresult <- nbinomTest(countdataset, condpair[1], condpair[2])
    filename <- paste0(fileprefix, condpair[1] , condpair[2])
    write.table( nbinomresult, file = filename, quote = FALSE, sep = '\t')
    return( list(condA = condpair[1], condB = condpair[2],
                 filename = filename))
}

f.deseq_testallpairs <- function(countdataset, fileprefix ) {
    FUN = function (x)  f.deseq_do_andwrite(countdataset,x,fileprefix)
    nbinomindex<- apply(X = combn(levels(conditions(countdataset)),m=2),
                        MARGIN = 2, FUN = FUN)
    return(nbinomindex)
}

load('./intermediate_results/cds.RData')
nbinomindex <- f.deseq_testallpairs(cds, './intermediate_results/nbinomresult')
save( list='nbinomindex', file='./intermediate_results/nbinomindex.RData')

