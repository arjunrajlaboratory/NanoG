dir.create(file.path('outputdata/figtophits'), showWarnings = FALSE)

library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)
library(stringr)
library(data.table)
library(grid)
library(scales)
Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')

load('./intermediate_results/m.deseq.RData')
load('./intermediate_results/m.rpkm.RData')
load('./intermediate_results/goannot.RData')


genelist <- list(
#     Ctrl = c('Actb'),
    Pluripotency = c('Myc','Stat3','Sox13','Sall4','Pou5f1','Sox2',
	'Klf4','Esrrb','Zfp42','Nanog','Pecam1','Fgf4'),
    Epi = c('Fgf5'),
    TE = c('Tfap2a','Cdx2','Cbx8','Gata3','Tead4'),
    PE = c('Gata4','Gata6','Pdgfra','Sox17','Fgfr2','Creb3l2'),
    Meso = c('T','Lefty1','Msx1','Tbx6'),
    Ectoderm = c('Nes','Pax3','Sox3','Crabp2')
)







# hitlist <- f.makehitlist( 'GFP.NEG', m.deseq, m.rpkm)
# hitspecial <- f.cross.wGenelist( hitlist, genelist)
# 
# yadjust <- numeric(length(hitspecial$GeneSymbol))
# names(yadjust) <- as.character(hitspecial$GeneSymbol)
# xadjust <- yadjust
# yadjust <- yadjust + 0.2 *sign(hitspecial$log2Fold )
# yadjust[c('Esrrb','Sox2')] <- 0.2
# yadjust[c('Pou5f1')] <- 0.2
# yadjust[c('Nanog')] <- -0.3
# yadjust[c('Msx1')] <- -0.2
# 
# hitspecial <- cbind(hitspecial, xadjust=xadjust, yadjust=yadjust)
# 
# q1 <- ggplot(hitlist, aes( x=log2(maxRPKM), y=log2Fold) ) + 
#     geom_point(size = 0.5) 
# q2 <- ggplot(hitspecial, aes(x=log2(maxRPKM), y=log2Fold)) + 
#     geom_point(size=2) +
#     geom_text( data = hitspecial, 
# 	aes(label = GeneSymbol, 
# 	    x = log2(maxRPKM) + xadjust,
# 	    y = log2Fold + yadjust),
# 	vjust = 0.5, color = 'black', size=3, fontface='bold' ) 



q <- f.plot_hitlist_genelist('GFP.NEG',genelist, m.deseq, m.rpkm)
q <- q + geom_abline(intercept=2.8, slope=-0.13) + guides(color='none')
qgfpneg <- q
q <- f.plot_hitlist_genelist('ALL.DIFF',genelist, m.deseq, m.rpkm)
q <- q + geom_abline(intercept=2.8, slope=-0.13) + guides(color='none')
qalldiff <- q
setEPS()
postscript('./outputdata/figtophits/tophits.fold.maxrpkm.eps',width=7,height=4)
pushViewport( viewport( layout=grid.layout(1,2)))
print(qgfpneg, vp = viewport( layout.pos.col = 1) )
print(qalldiff, vp = viewport( layout.pos.col = 2 ) )
dev.off()




tophits <- ldply( list('GFP.NEG','ALL.DIFF'), 
    function(x) f.makehitlist(x, m.deseq, m.rpkm) )
tophits <- transform( tophits, istophit = (log2Fold > 2.8-0.13*log2(maxRPKM)) )


goannotlong <- f.golist.to.godf( goannot )

GOcatlist <- list( 
    GO.development = c(
	'GO:0009790' # embryo development
	, 'GO:0048856' # anatomical structure development
	, 'GO:0007389'  # pattern specification process 
	),
    GO.meiosis = c('GO:0007126'),
    GO.FGF = c('GO:0008543'),
    GO.Wnt = c('GO:0016055')
)
for (n in seq(along=GOcatlist) ){
    tophits[[ names(GOcatlist)[n] ]] <- tophits$GeneSymbol %in%
	as.character( subset( goannotlong, Category %in% GOcatlist[[n]])$GeneSymbol )
}


fractionDevHitsInGFPNEGByHitType <- ddply(subset(tophits, meta == 'GFP.NEG'), 
    c('istophit','signedsig'), 
    summarize, fractionDevelopmentalHits = mean(GO.development))

capture.output(print(fractionDevHitsInGFPNEGByHitType),
    file = 'outputdata/figtophits/developmentFractions.txt', 
    append = FALSE)

q <- ggplot(data=NULL) + geom_abline(intercept = 2.8, slope = -0.13)
q <- q + geom_abline(intercept = 0, slope=0)
q <- q + ylim(-2.5,5) + xlim(-5,10)
q <- q + xlab('log2(maxRPKM)') + ylab('log2Fold')
qempty <- q



setEPS()
postscript('./outputdata/figtophits/empty.tophits.eps', width = 7/4, height=2)
print(qempty)
dev.off()

