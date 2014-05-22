dir.create(file.path('outputdata/figHox'), showWarnings = FALSE)

library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)
library(stringr)
library(data.table)
library(grid)
Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')



base_size <- 10



load('./intermediate_results/m.deseq.RData')



dftoplot <- subset( m.deseq, meta %in% c('GFP.NEG', 'ALL.DIFF'))
dftoplot <- dcast( dftoplot, 'GeneSymbol + meta ~ variable' )
dftoplot <- transform(dftoplot, meta = factor( meta, c('GFP.NEG','ALL.DIFF')))




f.drawhox <- function(){
    dfin <- subset(dftoplot, grepl('^Hox[ab]',GeneSymbol))
    qhoxab <- f.logfold.dotplot( dfin) + xlim(-3,5)
    dfin <- subset(dftoplot, grepl('^Hox[cd]',GeneSymbol))
    qhoxcd <- f.logfold.dotplot( dfin) + xlim(-3,5)
    pushViewport(viewport(layout=grid.layout(1,2)))
    pushViewport(viewport(layout.pos.col = 1))
    print(qhoxab, vp = current.viewport())
    upViewport()
    pushViewport(viewport(layout.pos.col = 2))
    print(qhoxcd, vp = current.viewport())
    upViewport()
}


setEPS()
postscript('./outputdata/figHox/dotplot.hox.eps',width=4,height=3.6)
f.drawhox()
dev.off()
