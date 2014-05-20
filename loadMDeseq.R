
if (!exists('m.deseq')){

    Stangle('s_functions.Rnw');source('s_functions.R')

    f.meltdeseq <- function(condA,condB,nbinomindex){
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
    load('./intermediate_results/nbinomindex.RData')
    m.deseq <- f.meltdeseq('GFP','NEG',nbinomindex)
    m.deseq <- rbind(m.deseq, f.meltdeseq('ALL','DIFF',nbinomindex))
    m.deseq <- rbind(m.deseq, f.meltdeseq('NEG','DIFF',nbinomindex))
    m.deseq <- rbind(m.deseq, f.meltdeseq('GFP','DIFF',nbinomindex))
    m.deseq <- rbind(m.deseq, f.meltdeseq('PLUSS','GFP',nbinomindex))

}