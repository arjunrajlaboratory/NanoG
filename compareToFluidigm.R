
library(plyr)
library(reshape2)
library(ggplot2)
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')

fluidigm <- read.csv('./inputdata/fluidigmAvg.csv')
fluidigm <- rename(fluidigm, c(row.names='lowerGeneSymbol'))

load('intermediate_results/m.deseq.RData')

rnaseq <- subset(m.deseq, meta == 'GFP.NEG' & variable %in% c('log2Fold'))
rnaseq <- dcast(rnaseq, 'GeneSymbol ~ variable')
rnaseq$lowerGeneSymbol <- tolower(rnaseq$GeneSymbol)

comparison <- merge(fluidigm, rnaseq, by='lowerGeneSymbol')

source('makeFluidigmGeneList.R')

comparison <- f.cross.wGenelist( comparison, genelist) 

base_size <- 10 

q <- ggplot(comparison, aes(x=fold_2iL, y=log2Fold))
q <- q + geom_point(aes(shape=L1))
q <- q + scale_shape(guide=FALSE)
q <- q + xlab('fluidigm 2i logFold')
q <- q + ylab('rnaseq 2i logFold')
q <- q + geom_text(aes(label=GeneSymbol),size=3, hjust=0, vjust=0)
q <- q + theme_minimal(base_size=base_size)

setEPS()
postscript('./outputdata/figFluidigm/vsfluidigm.eps',width=6,height=6)
print(q)
dev.off()

