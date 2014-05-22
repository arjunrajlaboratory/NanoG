<<echo=FALSE>>=
options(continue="  ")
options(width=70)
@

\section{Normalization data for RPKM}

This chunk gets the read totals that should be used to normalize the
counts. These are found at the top of the \verb+feature_quant...+ files,
and I extracted them with shell scripts

<<>>=
library(stringr)
readtotals <- read.delim( "./inputdata/readsummaryfeaturequant", sep = ':',
                         header = FALSE)
readtotals[,1]<- str_sub(readtotals[,1], start=1, end=2)
readtotalsforRPKM <- readtotals[,2]
names(readtotalsforRPKM) <- str_sub(readtotals[,1], start=1, end=2)
save( list = 'readtotalsforRPKM',
     file = './intermediate_results/readtotalsforRPKM.RData')
@
