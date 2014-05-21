

library(plyr)
library(ggplot2)
library(reshape2)
library(stringr)

Stangle('s_functions.Rnw');source('s_functions.R')

indexMorgani <- read.delim('./externaldata/Morgani13/indexGSE45182.txt', 
    header=FALSE, stringsAsFactors=FALSE)
filedir <- './externaldata/Morgani13/'

names(indexMorgani) <- c('gsm', 'media', 'HV', 'none', 'rep')
indexMorgani <- within( indexMorgani, {
    rep <- paste0('rep', rep)
    HV <- ifelse(HV == 'HV-', 'HVneg', HV)
    HV <- ifelse(HV == 'HV+', 'HVpos', HV)
    media <- str_replace(media, '/', '')
    })
indexMorgani[['none']] <- NULL

datfiles <- list.files( filedir, pattern = '^GSM')
datMorgan <- list()
for (filename in datfiles){
    gsm <- str_extract(filename, '^GSM[[:digit:]]*')
    datMorgan[[gsm]] <- read.delim( paste0(filedir,filename) )
    names(datMorgan[[gsm]]) <- c('Refseq', 'rpkm')
}
datMorgan <- melt(datMorgan, id.vars=c('Refseq','rpkm'))
names(datMorgan) <- c('Refseq', 'rpkm', 'gsm')


# Load the Refseq annotation data in the same GSE series

refdata <- read.delim(paste0( filedir, 'GSE45182_RefSeqInformation.txt'),
    header=TRUE)
refdata[[4]] <- NULL
names(refdata) <- c('Refseq', 'length', 'GeneSymbol')
datMorgan <- merge( datMorgan, refdata, by = 'Refseq')


# Turns out they have already chosen only one isoform per gene and they
# are not repeated in their data tables

print(any(duplicated( subset(datMorgan, gsm=='GSM1098624', select='GeneSymbol'))))



datMorgan <- merge(datMorgan, indexMorgani, by ='gsm')


# According to the way that Morgani et al. normalized their RPKM, it makes
# the most sense to compute the average of the log2Folds obtained for each
# replicate between the two conditions.


m.morgani <- melt(datMorgan, measure.vars='rpkm')
m.morgani <- dcast(m.morgani, 'GeneSymbol + media + rep ~ HV + variable')
m.morgani <- within(m.morgani, HVneg.HVpos_log2Fold <- log2(HVpos_rpkm/HVneg_rpkm))
m.morgani <- melt(m.morgani, id.vars = c('GeneSymbol','media','rep'))
m.morgani <- dcast(m.morgani, 'GeneSymbol + media + variable ~ rep')
m.morgani <- within(m.morgani, {avg <- (rep1 + rep2)/2; diff <- (rep2-rep1)})
m.morgani <- melt(m.morgani, c('GeneSymbol','media','variable'), variable.name='rep')
m.morgani <- subset(m.morgani, rep %in% c('avg','diff'))
m.morgani <- f.colsplitandmerge( m.morgani, 'variable', '_', c('HV','variable'))
m.morgani <- within(m.morgani, meta <- paste(HV,media,rep,sep='_'))
m.morgani <- subset(m.morgani, !is.na(value))
m.morgani <- m.morgani[,c('meta','GeneSymbol','variable','value')]
m.morgani <- cbind(source = 'Morgani', m.morgani)



save(file = './intermediate_results/m.morgani.RData', list='m.morgani')
