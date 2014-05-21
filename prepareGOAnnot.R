

library(goseq)

load('./intermediate_results/fullTable.RData')

goannot <- getgo( as.character(fullTable$Gene_Symbol), 'mm9', 'geneSymbol',
    fetch.cats = c('GO:BP'))
save( file = './intermediate_results/goannot.RData', list = 'goannot')

