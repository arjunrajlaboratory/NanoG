
library(plyr)
library(reshape2)
library(stringr)
library(DESeq)
Stangle('s_functions.Rnw'); source('s_functions.R')



# \section{Converting our data to an RPKM basis}

# This requires that we first have to convert our data to RPKM
# - reads per kilobase per million mapped reads.


load('./intermediate_results/fullTable.RData')
load('./intermediate_results/cds.RData')

# \section{Normalization data for RPKM}

# This chunk gets the read totals that should be used to normalize the
# counts. These are found at the top of the \verb+feature_quant...+ files,
# and I extracted them with shell scripts


readtotals <- read.delim( "./inputdata/readsummaryfeaturequant", sep = ':',
                         header = FALSE)
readtotals[,1]<- str_sub(readtotals[,1], start=1, end=2)
readtotalsforRPKM <- readtotals[,2]
names(readtotalsforRPKM) <- str_sub(readtotals[,1], start=1, end=2)





cds_rpm <- newCountDataSet(countData=counts(cds),conditions(cds))
# Normalize by a million reads
sizeFactors(cds_rpm) <- readtotalsforRPKM / 10^6
cds_rpm <- estimateDispersions(cds_rpm, method='per-condition',
                                sharingMode='maximum',fitType='parametric')



# We want a function that for a given gene gives the estimated RPKM and
# the estimated error in the RPKM. In our case we have a few
# measurements for each condition, each with different library sizes
# $s_i$. The mean rpkm for that condition is $\hat{q}=1/m \sum k_i/s_i$
# where $k_i$ are the counts for that gene, and m is the total number of
# replicates for that condition.

# On the other hand, the variances for each read are given by $\sigma^2_{k_i} = s_i q
# + \alpha s_i^2 q^2$. Therefore, the variance on $\hat{q}$ is given by:
# \begin{eqnarray}
# \hat{\sigma}^2_{\hat{q}} &=& \frac{1}{m^2} \sum\sigma^2_{k_i}/s_i^2
# \nonumber \\
# &=& \frac{1}{m^2} \sum \frac{s_i\hat{q} + \alpha s_i^2\hat{q}^2
# }{s_i^2} \nonumber \\
# &=& \frac{1}{m^2} \sum \frac{\hat{q}}{s_i} + \alpha \hat {q}^2 \nonumber
# \end{eqnarray}

# This is implemented in the function:

# f.getrpkmstats



# Implement this below and save the RPKM stats for our dataset in a
# file.

ourRPKMvals <- lapply(dispTable(cds_rpm),
                      FUN = function (x) f.getrpkmstats(cds_rpm, x,
                      fullTable$Length))
ourRPKMvals <- do.call('cbind',ourRPKMvals)




m.rpkm <- ourRPKMvals
m.rpkm$GeneSymbol <- row.names(m.rpkm)
m.rpkm <- melt(m.rpkm, id = 'GeneSymbol')
m.rpkm <- cbind( m.rpkm, 
    f.fixedcolsplit( m.rpkm$variable , '.', c('meta','var'))) 
m.rpkm$source <- 'deseq'
m.rpkm <- m.rpkm[,c('source','meta','GeneSymbol','var','value')]
m.rpkm <- rename( m.rpkm, c( var = 'variable' ) )
save(file = './intermediate_results/m.rpkm.RData', list = 'm.rpkm')