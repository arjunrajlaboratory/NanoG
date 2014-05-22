dir.create(file.path('outputdata/figCoh'), showWarnings = FALSE)

library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)
library(stringr)
library(data.table)
library(grid)
library(gridExtra)
Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')



base_size <- 10 
load('./intermediate_results/goannot.RData')
goannotlong <- f.golist.to.godf(goannot)



load('./intermediate_results/m.deseq.RData')
dftoplot <- m.deseq
dftoplot <- subset(dftoplot, variable %in% c('log2Fold','issig','signedsig'))
dftoplot <- subset(dftoplot, meta %in% c('GFP.NEG','ALL.DIFF'))
dftoplot <- dcast(dftoplot, 'GeneSymbol ~ meta + variable')
dftoplot[,-1] <- sapply(dftoplot[,-1], function(x) ifelse(is.finite(x),x,NA) )


f.scatterprettierCoh <- function(gplotin, countlim = 60){
    gplotout <- gplotin + coord_equal()
#     gplotout <- gplotout +
# 	geom_rect(data=NULL,xmin=0,ymin=0,ymax=Inf,xmax=Inf, fill='#E69F00',alpha=0.1)+
# 	geom_rect(data=NULL,xmin=-Inf,ymin=0,ymax=Inf,xmax=0, fill='#CC79A7',alpha=0.1)+
# 	geom_rect(data=NULL,xmin=0,ymin=-Inf,ymax=0,xmax=Inf, fill='#CC79A7',alpha=0.1)+
# 	geom_rect(data=NULL,xmin=-Inf,ymin=-Inf,ymax=0,xmax=0, fill='#E69F00',alpha=0.1)
    gplotout <- gplotout + 
# 	geom_hline(yintercept = 0 ) + 
# 	geom_vline(xintercept = 0) +
	geom_bin2d(binwidth=c(0.2,0.2)) 
    rain.colors <- rev( rainbow(7)[1:6])
    gplotout <- gplotout + scale_fill_gradientn(
	colours=rain.colors,
	breaks = c(0,countlim),
	guide = guide_colourbar(direction='horizontal',title=NULL,
	    label.position='bottom', barwidth = 2, barheight=0.25,
	    ticks = FALSE, )
    )
    gplotout <- gplotout + theme_grey(base_size=base_size)
    gplotout <- gplotout + theme(
	legend.position = c(0.80,0.95),
	legend.background = element_blank(),
	plot.margin=unit(rep(0,4),'lines'),
	axis.text = element_text(size=rel(0.8), colour='black'),
	axis.ticks = element_line(colour='black'),
	axis.title.x = element_blank(),
	axis.title.y = element_blank()
# 	panel.grid.major = element_blank()
    )
    gplotout <- gplotout + scale_y_continuous(limits=c(-4,4))
    gplotout <- gplotout + scale_x_continuous(limits=c(-3,4))
    return(gplotout)
}
df2 <- subset( dftoplot, ALL.DIFF_issig & GFP.NEG_issig ) 
q <- ggplot( df2, aes(GFP.NEG_log2Fold, ALL.DIFF_log2Fold)) 
q <- f.scatterprettierCoh(q,40)
qcohscatter <- q



tabletoplot <- summarize( df2, 
    ALL.DIFF_highin = factor(ALL.DIFF_signedsig, levels=c(1,-1),
	labels=c('Diff','Stem')),
    GFP.NEG_highin = factor(GFP.NEG_signedsig, levels=c(-1,1), 
	labels=c('VNP+','VNP-')))
tabletoplot <- as.data.frame(table(tabletoplot))
tabletoplot <- dcast( tabletoplot, 'ALL.DIFF_highin ~ GFP.NEG_highin', 
    value.var='Freq')
row.names(tabletoplot) <- tabletoplot$ALL.DIFF_highin
tabletoplot[,1] <- NULL
cohtable <- tableGrob( tabletoplot, 
    gpar.coltext = gpar(cex=1,fontsize=base_size*0.8,fontface='bold'),
    gpar.rowtext = gpar(cex=1,fontsize=base_size*0.8,fontface='bold'))



f.drawcoh <- function(){
    vpcohscatter <- viewport(x=0.25,y=0.5,width=0.5,height=1,name='vpcohscatter')
    vpcohtable <- viewport(x = 0.75, y=0.25,width=0.5,height=0.5,name='vpcohtable')
    print(qcohscatter, vp=vpcohscatter)
    pushViewport(vpcohtable)
    grid.draw(cohtable)
    upViewport()
}



setEPS()
postscript('./outputdata/figCoh/Cohscatter.eps',width=4,height=2.4)
print(qcohscatter)
dev.off()





require(GO.db)
df2 <- transform( df2, 
    ALL.DIFF = factor(ALL.DIFF_signedsig, levels=c(1,-1),
	labels=c('Diff','Stem')),
    GFP.NEG = factor(GFP.NEG_signedsig, levels=c(-1,1), 
	labels=c('VNP+','VNP-')))
CohGOcats <- ddply(df2, c('GFP.NEG','ALL.DIFF'),
    function(x) f.makegreedyGOtable( genehits = x$GeneSymbol, 
	gotable = goannotlong, gosigdata = NULL, gosizecutoff = 750 ),
	.progress = 'text')
CohGOcats.top <- ddply(CohGOcats, c('GFP.NEG','ALL.DIFF'),
    head, n=5)


setEPS()
postscript('./outputdata/figCoh/CohGO_750.eps',width=7,height=7)
pushViewport(viewport(gp=gpar(fontsize=8)))
grid.draw(tableGrob(CohGOcats.top[,c(1:7)],show.rownames=FALSE))
dev.off()


CohGOcats <- ddply(df2, c('GFP.NEG','ALL.DIFF'),
    function(x) f.makegreedyGOtable( genehits = x$GeneSymbol, 
	gotable = goannotlong, gosigdata = NULL, gosizecutoff = 250 ),
	.progress = 'text')
CohGOcats.top <- ddply(CohGOcats, c('GFP.NEG','ALL.DIFF'),
    head, n=5)
setEPS()
postscript('./outputdata/figCoh/CohGO_250.eps',width=7,height=7)
pushViewport(viewport(gp=gpar(fontsize=8)))
grid.draw(tableGrob(CohGOcats.top[,c(1:7)],show.rownames=FALSE))
dev.off()



