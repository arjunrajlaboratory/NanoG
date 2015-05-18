

library(reshape2)
library(goseq)
library(plyr)

load('./intermediate_results/goannot.RData')
load('./intermediate_results/m.deseq.RData')


f.rungoseq <- function( geneids, member.of.subset, dataforbias, 
    gene2cat, test.cats='GO:BP', method = 'Wallenius'){
    DEgenes <- as.vector(member.of.subset)
    names(DEgenes) <- geneids
    bias.data <- as.vector(dataforbias)
    names(bias.data) <- geneids
    pwf <- nullp(DEgenes, 'mm9', 'geneSymbol') #, bias.data=bias.data)
    GOdat <- goseq(pwf, genome = 'mm9', id = 'geneSymbol',
    gene2cat = gene2cat, test.cats = 'GO:BP', method = method)
    GOdat$GOTERM <- Term(GOdat$category)
    return(GOdat)
}

require(GO.db)
dattoGO <- subset( m.deseq, meta %in% c('GFP.NEG','ALL.DIFF')) 
dattoGO <- dcast( dattoGO, 'GeneSymbol ~ meta + variable' )
genehits <- with(dattoGO, as.character(
    GeneSymbol[ ALL.DIFF_issig == 1 | GFP.NEG_issig == 1] ))
dattoGO$avgbaseMean <- with(dattoGO, (GFP.NEG_baseMean + ALL.DIFF_baseMean)/2 )
dattoGO <- subset( dattoGO, avgbaseMean > 0 )

f.calccats_hits <- function(){
    cats_hits <- with(dattoGO, 
    f.rungoseq(GeneSymbol, GFP.NEG_issig == 1 | ALL.DIFF_issig == 1, 
        log10(avgbaseMean), goannot, test.cats='GO:BP'))
    cats_hits <- rename(cats_hits, c(category='Category')) 
    return(cats_hits)
}

cats_hits <- f.calccats_hits()

save(list='cats_hits',file='./intermediate_results/cats_hits_new.RData')
