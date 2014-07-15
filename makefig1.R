dir.create(file.path('outputdata/fig1'), showWarnings = FALSE)

library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)
library(stringr)
library(data.table)
library(grid)
Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')

load('intermediate_results/m.deseq.RData')

base_size <- 10


sortdata <- list()
for (sortname in c('A','C')){
    sortdata[[sortname]] <- read.delim(paste0('./inputdata/sorter',sortname,'.txt'),header=TRUE)
    sortdata[[sortname]][['FITC']] <- sortdata[[sortname]][['FITC.E14']]
}


sortdata <- melt( sortdata, c('FITC'))
sortdata <- f.colsplitandmerge( sortdata, 'variable', '.', c('variable','meta'))
sortdata <- subset(sortdata, variable %in% c('Events','Fluorescence'))

dflabels <- data.frame( meta = c('NEG','DIFF','GFP','ALL','PLUSS') )
dflabels <- transform( dflabels, label = f.sample.labeller(value = meta))
dflabels$x <- c(3.2,3.2,2.2,3.2,2.5)
dflabels$y <- rep(2.5,5)

dftoplot <- subset( sortdata, L1 == 'A' )
dftoplot <- dcast(dftoplot, '... ~ variable')
dftoplot <- subset( dftoplot, Fluorescence > 0 )
dftoplot <- ddply(dftoplot, c('meta','L1'), transform, Freq=Events/sum(Events))

dfbackg <- subset(dftoplot, meta == 'E14', 
    select=c('FITC','Events','Fluorescence','Freq'))
dftoplot <- subset(dftoplot, meta != 'E14')
dftoplot$meta <- factor(dftoplot$meta, dflabels$meta)
dfstem.backg <- subset(dftoplot, meta == 'ALL')
dfstem.backg$meta <- 'DIFF'



q <- ggplot(data = dftoplot, aes(x = log10(Fluorescence), weight = Freq)) + 
    geom_density( data = dfbackg, fill='gray',color='gray', adjust=1/10) +
    geom_density( data = dfstem.backg, fill = 'white', color = 'gray', adjust = 1/10)+
    geom_density(adjust=1/10, fill='black') +
    geom_text( data = dflabels, aes(x=x,y=y,label = label, weight = NULL),size=4) +
    facet_wrap('meta',ncol=2, as.table=FALSE ) 
q <- q + scale_x_continuous(breaks=c(2,3,4),limits=c(1.5,4)) 
q <- q + scale_fill_discrete(guide = 'none')
q <- q + theme_minimal(base_size = base_size)
q <- q + theme(strip.text.x = element_blank(), strip.background = element_blank())


setEPS()
postscript('./outputdata/fig1/sorterA.eps',width=2.6,height=3.5)
print(q)
invisible(dev.off())



source('makeFluidigmGeneList.R')


dftoplot <- rbind( m.deseq )
dftoplot <- subset( dftoplot, variable %in% c('log2Fold','issig') )
dftoplot <- dcast( dftoplot, 'source+meta+GeneSymbol ~ variable')
dftoplot <- f.cross.wGenelist( dftoplot, genelist) 
df2 <- subset(dftoplot, meta %in% c('GFP.NEG','ALL.DIFF'))
df2$meta <- factor(df2$meta, c('GFP.NEG','ALL.DIFF'))



q <- ggplot(df2, aes(log2Fold, GeneSymbol, fill=meta )) +
#     geom_point(data=subset(df2,issig=='ishit'),
# 	color='black', size=3,fill='white') +
    geom_point(shape=21) + geom_vline(xintercept = 0) +
    geom_point( data = subset(df2,issig==1), shape= '.', color='red') +
    facet_grid( L1 ~ . , scales = 'free_y', space = 'free_y') +
    scale_fill_manual(values=c('black','white')) 
q <- q + theme_minimal(base_size=base_size) 
q <- q + theme(
    legend.position='top',
    strip.background = element_blank(),
    panel.margin = unit(0.5,'lines')
)
q <- q + guides( fill = guide_legend( ncol = 1, label.hjust = 0, title = NULL))
q <- q + ylab(NULL) + xlab('log2Fold')


setEPS()
postscript('./outputdata/fig1/abranchesdotplot.eps',width=1.8,height=7) 
print(q)
invisible(dev.off())



