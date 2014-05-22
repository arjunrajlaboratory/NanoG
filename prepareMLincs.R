
library(stringr)
library(plyr)
library(reshape2)
Stangle('s_functions.Rnw');source('s_functions.R')

f.meltlincs <- function(condA,condB,nbinomindex){
    nb <- f.loadntest(condA,condB,nbinomindex)
    nb <- nb[,c('id','baseMean','log2FoldChange','pval','padj')]
    nb$issig <- f.meetsfdr( nb)
    nb$signedsig <- ifelse( nb$issig & nb$log2FoldChange > 0, 1, 0)
    nb$signedsig <- ifelse( nb$issig & nb$log2FoldChange < 0, -1, nb$signedsig)
    nb <- rename( nb, c( id = 'GeneSymbol',log2FoldChange = 'log2Fold' ) )
    nb <- melt( nb, id = 'GeneSymbol')
    nb <- cbind( source= 'deseq', meta = paste0(condA,'.',condB), nb)
    return(nb)
}


load('./intermediate_results/nbinomindex.wlincs.RData')
m.lincs <- f.meltlincs('GFP','NEG',nbinomindex.wlincs)
m.lincs <- rbind(m.lincs, f.meltlincs('ALL','DIFF',nbinomindex.wlincs))
m.lincs <- rbind(m.lincs, f.meltlincs('NEG','DIFF',nbinomindex.wlincs))
m.lincs <- rbind(m.lincs, f.meltlincs('GFP','DIFF',nbinomindex.wlincs))
m.lincs <- rbind(m.lincs, f.meltlincs('PLUSS','GFP',nbinomindex.wlincs))
save( file = './intermediate_results/m.lincs.RData', list = 'm.lincs')
