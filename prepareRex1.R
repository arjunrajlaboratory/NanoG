

library(plyr)
library(reshape2)

load('./intermediate_results/fullTable.RData')

datasetnames <- c('GSM850398','GSM850399')
results <- list()
for (dataset in datasetnames){
    x <- read.delim( paste0( './externaldata/marksetalcell/transcripts_', dataset),
     header = FALSE )
    names(x) <- c('transcript','chrpos','Ucounts','NUcounts','Length','strand','first_annot')
    x <- merge(x, fullTable[,c('first_annot','Gene_Symbol')], by = 'first_annot',
    all.x=FALSE, all.y=TRUE)
    totalreads <- read.delim( paste0('./externaldata/marksetalcell/features_',dataset),
    header = FALSE, nrows=1, sep = ' ')
    totalreads <- totalreads[1,3]
    x[['rpkm']] <- x$Ucounts / (x$Length/10^3) / (totalreads/10^6)
    x <- rename( x, c(Gene_Symbol = 'GeneSymbol', Ucounts = 'counts'))
    results[[dataset]] <- x[,c('GeneSymbol','counts','rpkm')]
}
results <- rename(results, c(GSM850398 = 'Rexneg', GSM850399 = 'Rexpos'))
m.rex1 <- melt(results, 'GeneSymbol') 
m.rex1 <- rename(m.rex1, c(L1 = 'meta'))
temp <- dcast( subset(m.rex1, variable == 'rpkm'), 'GeneSymbol ~ meta')
temp$log2Fold <- log2(temp$Rexneg / temp$Rexpos)
temp <- melt( temp[,c('GeneSymbol','log2Fold')], 'GeneSymbol')
m.rex1 <- rbind( m.rex1, cbind(meta = 'Rexpos.Rexneg', temp))
m.rex1 <- cbind(source='marks', m.rex1)
m.rex1 <- m.rex1[,c('source','meta','GeneSymbol','variable','value')]
save(file='./intermediate_results/m.rex1.RData', list = 'm.rex1')


