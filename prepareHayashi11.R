# \section{Import data from Hayashi et al. cell 2011}



library(plyr)
library(reshape2)
library(stringr)

gseraw <- read.delim('./externaldata/HayashiCell11/GSE30056.txt.labelled',header=FALSE)
names(gseraw) <- c('gsm','ID','intensity','qc')
# log-transform the intensity:
gseraw <- transform( gseraw, intensity = log2(intensity))
gpl <- read.delim('./externaldata/HayashiCell11/GPL1261-14790.txt',skip=16,header=TRUE)
gse <- merge(gseraw, gpl[,c('ID','Gene.Symbol')], by='ID')
gse <- rename(gse, c(Gene.Symbol='GeneSymbol'))
gse <- subset(gse, GeneSymbol != '' )
table( str_detect(gse$GeneSymbol, '///'))
gse <- subset(gse, !str_detect(GeneSymbol, '///'))



sample.id <- read.delim('./externaldata/HayashiCell11/indexGSE30056.txt',header=FALSE)
sample.id <- sample.id[,c(1:3)]
names(sample.id) <- c('gsm', 'cond','rep') 
gse <- merge(gse, sample.id, by='gsm')


# For each probe, average across biological replicates. 
# Then, for each gene symbol, average over all probes that
# correspond to it, if there are more than 1. 
# The variance in the first step is calculated for each probe

library(data.table)
d.gse <- data.table(gse)
d.gse <- d.gse[, list( avgexpress = mean(intensity), 
    varexpress = var(intensity)), by='cond,ID,GeneSymbol']
d.gse <- d.gse[, list( avgexpress = mean(avgexpress),
    varexpress = mean(varexpress)), by='cond,GeneSymbol']



d.gse <- melt(d.gse, id=c('cond','GeneSymbol'))
d.gse <- rename( d.gse, list(cond = 'meta'))
d.gse$source <- 'Hayashi11'
d.gse <- d.gse[,c('source','meta','GeneSymbol','variable','value')]




# function that finds the log2Fold changes to ESC 
f.baseline <- function(df) {
    df <- subset(df, df$variable == 'avgexpress' )
    df <- dcast(df, 'source+GeneSymbol ~ meta')
    df$baseline <- df$ESC
    df <- melt(df, id = c('source','GeneSymbol','baseline'))
    df <- transform( df, value = value - baseline)
    df <- rename(df, c(variable = 'meta'))
    df$meta <- paste0('ESC.',df$meta)
    df$variable <- 'log2Fold'
    df[,c('source','meta','GeneSymbol','variable','value')]
}



m.hayashi <- rbind(d.gse, f.baseline(d.gse))
save(file = './intermediate_results/m.hayashi.RData', list ='m.hayashi')

