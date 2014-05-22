

Stangle('s_functions.Rnw');source('s_functions.R')
library(reshape2)
library(plyr)
library(stringr)


# \section{Import Microarray data}
# Read in the SOFT-formatted data for GSE40355, after passing
# through a perl script that adds a column label indicating which
# Sample (GEO accession code) each comes from.

# The sample annotations are from the GSE40335 series GEO page,
# labelled by me to have names consistent with those used
# throughout their paper. 

# The gene annotations corresponding to the microarray probes
# come from the GPL6246 platform data set at GEO. 

d.mac <- read.delim('./externaldata/MacArthur12/GSE40335.txt.labelled',header=FALSE)
names(d.mac) <- c('gsm','probeID','intensity')
head(d.mac,n=3)
mac.names <- read.delim('./externaldata/MacArthur12/sampledescriptions.txt', header=FALSE)
names(mac.names) <- c('gsm','cond','fullspec')
mac.names$rep <- str_extract(mac.names$fullspec, 'Rep.')
mac.names$cond <- paste0('t',mac.names$cond)
head(mac.names,n=3)
d.annot <- read.delim('./externaldata/MacArthur12/GPL6246.annot.IDandGeneSymbol', header=FALSE)
names(d.annot) <- c('probeID','GeneSymbol')
head(d.annot,n=3)


# Remove probes that do not match any gene symbol.
# Also remove probes labelled for multiple gene symbols
# (denoted by a verb+///+ separator)

d.annot <- subset(d.annot, GeneSymbol != '' )
table( str_detect(d.annot$GeneSymbol, '///'))
d.annot <- subset(d.annot, !str_detect(GeneSymbol, '///'))
d.mac <- merge(d.mac, d.annot, by='probeID', all.x=TRUE, all.y =FALSE)
d.mac <- subset(d.mac, !is.na(GeneSymbol))
d.mac <- cbind( d.mac,  mac.names[ match( d.mac$gsm, mac.names$gsm ),c('cond','rep')])
head(d.mac, n=3)


# For each probe, average across biological replicates. 
# Then, for each gene symbol, average over all probes that
# correspond to it, if there are more than 1. 
# The variance in the first step is calculated for each probe

table(table(as.character(d.annot$GeneSymbol)))
library(data.table)
d.mac <- data.table(d.mac)
m.mac <- d.mac[, list( avgexpress = mean(intensity), 
    varexpress = var(intensity)), by='cond,probeID,GeneSymbol']
m.mac <- m.mac[, list( avgexpress = mean(avgexpress),
    varexpress = mean(varexpress)), by='cond,GeneSymbol']


# Melt the data into our common format.

m.mac <- melt(m.mac, id=c('cond','GeneSymbol'))
m.mac <- rename( m.mac, list(cond = 'meta'))
m.mac$source <- 'macarthur'
# Order the factors in the same order as in their paper
m.mac$meta <- factor(m.mac$meta, c('t0h','t24h','t36h','t36h.R','t72h',
    't84h','t84h.R','t120h','t132h','t132h.R'))
m.mac <- m.mac[,c('source','meta','GeneSymbol','variable','value')]


# Before saving, add a baseline-subtracted version of the data, giving
# the log2Fold of each gene from the \verb+t0h+ data.


# function that finds the log2Fold changes to t0h 
f.macbaseline <- function(df) {
	df <- subset(df, df$variable == 'avgexpress' )
	df <- dcast(df, 'source+GeneSymbol ~ meta')
	df$baseline <- df$t0h
	df <- melt(df, id = c('source','GeneSymbol','baseline'))
	df <- transform( df, value = value - baseline)
	df <- rename(df, c(variable = 'meta'))
	df$variable <- 'log2Fold'
	df[,c('source','meta','GeneSymbol','variable','value')]
}
m.mac <- rbind( m.mac, f.macbaseline(m.mac))


save(file = './intermediate_results/m.mac.RData', list = 'm.mac')

