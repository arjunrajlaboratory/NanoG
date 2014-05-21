
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(grid)
Stangle('s_functions.Rnw');source('s_functions.R')

base_size <- 10 


load('./intermediate_results/m.deseq.RData')


dftoplot <- dcast(m.deseq, 'GeneSymbol+meta ~ variable')
dftoplot$issig <- factor( dftoplot$issig, c(1,0))
dftoplot <- subset(dftoplot, meta %in% c('GFP.NEG','PLUSS.GFP'))
dftoplot$meta <- factor( dftoplot$meta, c('GFP.NEG','PLUSS.GFP'), 
    labels = c(' \n ','\n '))
q <- ggplot(dftoplot, aes(x=log2Fold, fill=issig))
q <- q + geom_histogram( binwidth = 0.1)
q <- q + facet_grid(.~meta, scales='free_y') 
q <- q + scale_fill_manual( values = c('black','grey90'), guide='none')
q <- q + xlim(-2,3)
q <- q + ylab('genes/bin\n')
q <- q + xlab('\n\n ')
q <- q + theme_minimal(base_size = base_size)
q <- q + theme(
    strip.background = element_blank(),
    plot.margin = unit(rep(0,4),'lines'),
    axis.title.x=element_text(size=rel(0.8)),
    axis.title.y=element_text(size=rel(0.8))
)



setEPS()
postscript('./outputdata/figPLUSS/foldPLUSS.eps',width=3,height=2)
print(q)
dev.off()

