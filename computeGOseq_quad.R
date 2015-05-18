

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

# Up in NEG, DIFF.
require(GO.db)
dattoGO <- subset( m.deseq, meta %in% c('GFP.NEG','ALL.DIFF')) 
dattoGO <- dcast( dattoGO, 'GeneSymbol ~ meta + variable' )
genehits <- with(dattoGO, as.character(
    GeneSymbol[ ALL.DIFF_signedsig == 1 | GFP.NEG_signedsig == 1] ))
dattoGO$avgbaseMean <- with(dattoGO, (GFP.NEG_baseMean + ALL.DIFF_baseMean)/2 )
dattoGO <- subset( dattoGO, avgbaseMean > 0 )

f.calccats_hits <- function(){
    cats_hits <- with(dattoGO, 
    f.rungoseq(GeneSymbol, GFP.NEG_signedsig == 1 | ALL.DIFF_signedsig == 1, 
        log10(avgbaseMean), goannot, test.cats='GO:BP'))
    cats_hits <- rename(cats_hits, c(category='Category')) 
    return(cats_hits)
}

cats_hits_UpNeg_UpDiff <- f.calccats_hits()

save(list='cats_hits_UpNeg_UpDiff',file='./intermediate_results/cats_hits_UpNeg_UpDiff.RData')

# Down in NEG, DIFF.
require(GO.db)
dattoGO <- subset( m.deseq, meta %in% c('GFP.NEG','ALL.DIFF')) 
dattoGO <- dcast( dattoGO, 'GeneSymbol ~ meta + variable' )
genehits <- with(dattoGO, as.character(
  GeneSymbol[ ALL.DIFF_signedsig == -1 | GFP.NEG_signedsig == -1] ))
dattoGO$avgbaseMean <- with(dattoGO, (GFP.NEG_baseMean + ALL.DIFF_baseMean)/2 )
dattoGO <- subset( dattoGO, avgbaseMean > 0 )

f.calccats_hits <- function(){
  cats_hits <- with(dattoGO, 
                    f.rungoseq(GeneSymbol, GFP.NEG_signedsig == -1 | ALL.DIFF_signedsig == -1, 
                               log10(avgbaseMean), goannot, test.cats='GO:BP'))
  cats_hits <- rename(cats_hits, c(category='Category')) 
  return(cats_hits)
}

cats_hits_DownNeg_DownDiff <- f.calccats_hits()

save(list='cats_hits_DownNeg_DownDiff',file='./intermediate_results/cats_hits_DownNeg_DownDiff.RData')



# Down in NEG, Up in DIFF.
require(GO.db)
dattoGO <- subset( m.deseq, meta %in% c('GFP.NEG','ALL.DIFF')) 
dattoGO <- dcast( dattoGO, 'GeneSymbol ~ meta + variable' )
genehits <- with(dattoGO, as.character(
  GeneSymbol[ ALL.DIFF_signedsig == 1 | GFP.NEG_signedsig == -1] ))
dattoGO$avgbaseMean <- with(dattoGO, (GFP.NEG_baseMean + ALL.DIFF_baseMean)/2 )
dattoGO <- subset( dattoGO, avgbaseMean > 0 )

f.calccats_hits <- function(){
  cats_hits <- with(dattoGO, 
                    f.rungoseq(GeneSymbol, GFP.NEG_signedsig == -1 | ALL.DIFF_signedsig == 1, 
                               log10(avgbaseMean), goannot, test.cats='GO:BP'))
  cats_hits <- rename(cats_hits, c(category='Category')) 
  return(cats_hits)
}

cats_hits_DownNeg_UpDiff <- f.calccats_hits()

save(list='cats_hits_DownNeg_UpDiff',file='./intermediate_results/cats_hits_DownNeg_UpDiff.RData')




# Up in NEG, Down in DIFF.
require(GO.db)
dattoGO <- subset( m.deseq, meta %in% c('GFP.NEG','ALL.DIFF')) 
dattoGO <- dcast( dattoGO, 'GeneSymbol ~ meta + variable' )
genehits <- with(dattoGO, as.character(
  GeneSymbol[ ALL.DIFF_signedsig == -1 | GFP.NEG_signedsig == 1] ))
dattoGO$avgbaseMean <- with(dattoGO, (GFP.NEG_baseMean + ALL.DIFF_baseMean)/2 )
dattoGO <- subset( dattoGO, avgbaseMean > 0 )

f.calccats_hits <- function(){
  cats_hits <- with(dattoGO, 
                    f.rungoseq(GeneSymbol, GFP.NEG_signedsig == 1 | ALL.DIFF_signedsig == -1, 
                               log10(avgbaseMean), goannot, test.cats='GO:BP'))
  cats_hits <- rename(cats_hits, c(category='Category')) 
  return(cats_hits)
}

cats_hits_UpNeg_DownDiff <- f.calccats_hits()

save(list='cats_hits_UpNeg_DownDiff',file='./intermediate_results/cats_hits_UpNeg_DownDiff.RData')
