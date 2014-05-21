
library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)
library(stringr)
library(data.table)
library(grid)
library(gridExtra)
library(GO.db)
Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')

base_size <- 10 


load('./intermediate_results/m.deseq.RData')
load('./intermediate_results/goannot.RData')
load('./intermediate_results/cats_hits.RData')


dftoplot <- dcast(m.deseq, 'GeneSymbol+meta ~ variable')
dftoplot$issig <- factor( dftoplot$issig, c(1,0))
dftoplot <- subset(dftoplot, meta %in% c('GFP.NEG','ALL.DIFF'))
dftoplot$meta <- factor( dftoplot$meta, c('GFP.NEG','ALL.DIFF'), 
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
qfold <- q



f.draw.arrowedlegend.x <- function( negxlab, posxlab ){
    ypos <- unit(-1.5, 'strheight', 'A')
    grid.text( negxlab, x=0.4, y=ypos, just=c('right','center'))
    grid.text( posxlab, x=0.6, y=ypos, just=c('left','center'))
    grid.lines( x=c(0.42,0.58),y=ypos,gp=gpar(fill='black',lwd=1),
	arrow=arrow(angle=45,type='closed', 
	    length=unit(0.5,'strwidth','A'), ends='both'))
}
f.drawfoldchange <- function(gplot,vp){
    print(gplot,vp=vp )
    
    vpdepth <- downViewport('axis-b.5-4-5-4')
    f.draw.arrowedlegend.x('VNP+','VNP-')
    subtitleypos <- unit(-3.5, 'strheight','A')
    grid.text( 'log2( RNA fold change)',
	x=1, y=subtitleypos, just=c('center','center'))
    upViewport(vpdepth)
    
    vpdepth <- downViewport('axis-b.5-6-5-6')
    f.draw.arrowedlegend.x('Stem','Diff')
    upViewport(vpdepth)
    
    vpdepth <- downViewport('panel.4-4-4-4')
    pushViewport(viewport( clip='off'))
    grid.rect( x=0.7,y=0.7,
	width= stringHeight('h'), height= stringHeight('h'),
	just=c('right','top'), gp=gpar(fill='black'))
    grid.text( 'hits at\n10%\nFDR', x=0.75, y=0.7,
	just=c('left','top'), gp=gpar(lineheight=0.9))
    upViewport(1)
    upViewport(vpdepth)
}



setEPS()
postscript('./outputdata/fig2/all.foldchanges.eps',width=5.6,height=3.5)
print(q)
invisible(dev.off())



# \section{Broad GO results}




f.golist.to.godf <- function(goannot){ 
    gonames <- unlist(goannot,use.names=FALSE)
    genenames <- ldply(goannot, length)
    genenames <- rep( genenames[['.id']], genenames[['V1']])
    goannotlong <- data.frame( GeneSymbol = genenames, Category = gonames)
    return( data.table(goannotlong) )
}

goannotlong <- f.golist.to.godf( goannot )





dattoGO <- subset( m.deseq, meta %in% c('GFP.NEG','ALL.DIFF')) 
dattoGO <- dcast( dattoGO, 'GeneSymbol ~ meta + variable' )
genehits <- with(dattoGO, as.character(
    GeneSymbol[ ALL.DIFF_issig == 1 | GFP.NEG_issig == 1] ))


GOcats.bygreedy <- f.makegreedyGOtable( genehits = genehits, 
    gotable = goannotlong, gosigdata = cats_hits, gosizecutoff = 750 )


dftoplot <- f.deseqXgreedyGO( m.deseq, GOcats.bygreedy, goannotlong )
dftoplot <- transform( dftoplot, 
    highin = factor(highin, labels = paste( 'genes up\n in', 
    f.sample.labeller(value=levels(highin)))))
dftoplot <- transform(dftoplot, TERM = factor(TERM, 
	labels = ifelse( levels(TERM) =='regulation of RNA metabolic process',
	    'regulation of RNA metab. proc.', levels(TERM))))
gfpneglabel <- 'VNP- relative to VNP+\nDown        Up'
alldifflabel <- 'Diff relative to Stem\nDown        Up'
dftoplot <- transform(dftoplot, meta = factor(meta, 
    levels=c('GFP.NEG','ALL.DIFF'),
    labels=c(gfpneglabel, alldifflabel)))


q <- f.plotFracsxdeseq.twoside( dftoplot )
q <- q + scale_y_continuous(breaks=c(-200,-100,0,100,200), 
    labels=c(200,100,0,100,200),limits=c(-250,250)) 
q <- q + scale_fill_manual(values=c('black','grey80'))
q <- q + geom_hline(yintercept=0,color='white')
# q <- q + scale_fill_brewer(palette='Set1')
q <- q + theme_minimal(base_size=base_size) + theme(
    strip.background = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size=rel(0.8))
)

q <- q + ylab('\nNumber of hit genes falling*\n into GO category')
q <- q + xlab(NULL)
q <- q + theme(plot.margin=unit(c(0, 0, 0, 0), "lines"))
qgoxdeseq <- q



dftoplot <- GOcats.bygreedy
dftoplot$issig <- factor( dftoplot$issig, c(FALSE, TRUE ) ) 
dftoplot$stripplaceholder <- ' \n '
dftoplot <- transform(dftoplot, TERM = factor(TERM, 
	labels = ifelse( levels(TERM) =='regulation of RNA metabolic process',
	    'regulation of RNA metab. proc.', levels(TERM))))
q <- ggplot( data = subset( dftoplot, TERM != 'biological_process'), 
    aes( x = factor(TERM,rev(levels(TERM) ))) ) 
q <- q + geom_point( aes( y = Frac.hits, size = issig ))
q <- q + scale_size_manual( values=c(1,2), drop=FALSE, guide ='none')
q <- q + geom_hline( yintercept = dftoplot$Frac.hits.rndm[1] )
q <- q + coord_flip() + facet_grid('.~stripplaceholder')
q <- q + scale_y_continuous(breaks=c(0.3,0.4,0.5),labels=c('.3','.4','.5')) 
q <- q + expand_limits(y=0.5)
q <- q + theme_minimal(base_size=base_size)
q <- q + theme(
    strip.background = element_blank(),
    plot.margin=unit(c(0, 0, 0, 0), "lines"),
    axis.title.x=element_text(hjust=1,size=rel(0.8)),
    axis.ticks.y=element_blank()
)
q <- q + xlab('')
q <- q + ylab('\nFraction of annotated\n genes that are hits')
qgotable <- q


f.drawGOanalysis <- function(){
    pushViewport(viewport(layout = grid.layout(1,3,
	widths = unit.c(stringWidth('regulation of RNA metab. proc.'), 
	    unit(c(0.2, 0.8), c('null','null')))
	    )))
    print(qgoxdeseq, vp = viewport( layout.pos.col=3))
    print(qgotable, vp = viewport( layout.pos.col=1:2))
    upViewport()
}


setEPS()
postscript('./outputdata/fig2/GOprojected.eps',width=7,height=2.4)
f.drawGOanalysis()
invisible(dev.off())


# \section{TF binding plots}

# <<>>=
# load('./intermediate_results/m.kim.RData')
# @

# <<>>=
# dftoplot <- m.kim[,c('meta','value')]
# dftoplot <- subset( dftoplot, !(meta == 'birA'))
# names(dftoplot) <- c('TF','GeneSymbol')
# dftoplot <- merge(dftoplot, dftoplot, by='GeneSymbol')
# # dftoplot <- transform( dftoplot, TF.x = factor(TF.x, levels(m.kim$meta) ),
# #     TF.y = factor(TF.y, levels(m.kim$meta)))
# dftoplot <- transform( dftoplot, 
#     TF.x = factor(toupper(TF.x), toupper(levels(m.kim$meta)) ),
#     TF.y = factor(toupper(TF.y), toupper(levels(m.kim$meta))))
# df2 <- subset( m.deseq, meta %in% c('ALL.DIFF','GFP.NEG' ) )
# df2 <- subset( df2, variable %in% c('log2Fold','issig','signedsig') )
# df2 <- dcast( df2, 'GeneSymbol + meta ~ variable')
# df3 <- merge( dftoplot, df2, by='GeneSymbol', all.x=FALSE, all.y=FALSE)
# require(data.table)
# dftoplot <- as.data.table(df3)
# dftoplot <- dftoplot[, list( Num = length(GeneSymbol)), 
#     by = c('TF.x','TF.y','meta','signedsig')]
# dftoplot <- dftoplot[, list( Num = Num, signedsig = signedsig, 
#     Frac = Num/sum(Num) ), 
#     by = c('TF.x','TF.y','meta')]
# dftoplot <- subset( dftoplot, !(signedsig == 0 ))
# dftoplot$signedsig <- factor(dftoplot$signedsig, c(-1,1))

# dfuniv <- as.data.table(df2)
# dfuniv <- dfuniv[, list( Num = length(GeneSymbol)),
#     by = c('meta','signedsig')]
# dfuniv <- dfuniv[, list( Num = Num, signedsig = signedsig, Frac=Num/sum(Num)),
#     by = c('meta')]
# dfuniv <- subset(dfuniv, !(signedsig == 0 ))
# dfuniv$signedsig <- factor(dfuniv$signedsig, c(-1,1))

# dftoplot <- subset( dftoplot, as.character(TF.x) == as.character(TF.y) )
# dftoplot$highin <- f.getHighIn( dftoplot )
# dftoplot$highin <- factor( dftoplot$highin, c('GFP','NEG','ALL','DIFF') )
# dftoplot <- transform( dftoplot, 
#     highin = factor(highin, labels = paste( 'up in\n', 
#     f.sample.labeller(value=levels(highin)))))
# dfuniv$highin <- f.getHighIn( dfuniv )
# dfuniv$highin <- factor( dfuniv$highin, c('GFP','NEG','ALL','DIFF') )
# dfuniv <- transform( dfuniv, 
#     highin = factor(highin, labels = paste( 'up in\n', 
#     f.sample.labeller(value=levels(highin)))))
# q <- ggplot( dftoplot, aes( y = Frac, x = TF.x ) ) +
#     geom_hline( data=dfuniv, aes(yintercept = Frac )) +
#     geom_point() + facet_grid('. ~ highin') + coord_flip()
# q <- q + theme_minimal(base_size=base_size) 
# q <- q + theme(
#     panel.margin = unit(0.5,'lines'),
#     axis.title.x = element_text(hjust=0.5,size=rel(0.8)),
#     axis.title.y = element_text(size=rel(0.8))
# )
# # q <- q + theme(plot.margin=unit(rep(0,4),'lines'))
# q <- q + scale_y_continuous(limits=c(0,0.3), breaks=c(0,0.1,0.2,0.3), 
#     labels = c('0','.1','.2',''))
# q <- q + ylab('\nFraction of genes bound by\ntranscription factor that are hits')
# q <- q + xlab('\nTranscription factor\n')
# #q <- q + theme(axis.text.x = element_text(angle = -90,hjust=0, vjust=0.5))
# qtfxdeseq <- q
# @
# <<>>=
# setEPS()
# postscript('./outputdata/fig2/TFxdeseq.eps', width = 4.5, height=2.4)
# print(q)
# invisible(dev.off())
# @

# \section{Scatter plots}

# <<>>=
# load('./intermediate_results/m.deseq.RData')
# load('./intermediate_results/m.rex1.RData')
# load('./intermediate_results/m.marksrna.RData')
# load('./intermediate_results/m.hayashi.RData')
# load('./intermediate_results/m.morgani.RData')
# @
# <<>>=
# dftoplot <- rbind(m.rex1, m.deseq, m.marksrna, m.hayashi, m.morgani)
# dftoplot <- subset(dftoplot, variable %in% c('log2Fold','issig'))
# dftoplot <- subset(dftoplot, meta %in% c('GFP.NEG','ALL.DIFF','Rexpos.Rexneg',
#     'twoi.serum','ESC.EpiLC:day1','HVneg.HVpos_2iLIF_avg'))
# dftoplot <- dcast(dftoplot, 'GeneSymbol ~ meta + variable')
# dftoplot[,-1] <- sapply(dftoplot[,-1], function(x) ifelse(is.finite(x),x,NA) )
# @

# <<>>=
# f.scatterprettier <- function(gplotin, countlim = 60){
#     gplotout <- gplotin + coord_equal()
#     gplotout <- gplotout + 
# # 	geom_hline(yintercept = 0 ) + 
# # 	geom_vline(xintercept = 0) +
# 	geom_bin2d(binwidth=c(0.2,0.2)) 
#     rain.colors <- rev( rainbow(7)[1:6])
#     gplotout <- gplotout + scale_fill_gradientn(
# 	colours=rain.colors,
# 	breaks = c(0,countlim),
# 	guide = guide_colourbar(direction='horizontal',title=NULL,
# 	    label.position='top', barwidth = 2, barheight=0.25,
# 	    ticks = FALSE, )
#     )
#     gplotout <- gplotout + theme_grey(base_size=base_size)
#     gplotout <- gplotout + theme(
# 	legend.position = c(0.60,0.12),
# 	legend.background = element_blank(),
# 	plot.margin=unit(rep(0,4),'lines'),
# 	axis.text = element_text(size=rel(0.8), colour='black'),
# 	axis.ticks = element_line(colour='black'),
# 	axis.title.x = element_blank(),
# 	axis.title.y = element_blank()
# # 	panel.grid.major = element_blank()
#     )
#     gplotout <- gplotout + scale_y_continuous(limits=c(-5,5),breaks=c(-4,-2,0,2,4))
# #     gplotout <- gplotout + geom_hline(yintercept=0) + geom_vline(xintercept=0)
#     return(gplotout)
# }
# genestokeep <- subset( m.deseq, variable == 'issig' & value == 1 & meta == 'GFP.NEG',
#     select = 'GeneSymbol' , drop = TRUE)
# df2 <- subset( dftoplot, GeneSymbol %in% genestokeep) 
# q <- ggplot( df2, aes(GFP.NEG_log2Fold, Rexpos.Rexneg_log2Fold))  
# q <- f.scatterprettier(q)
# q <- q + xlim(-2,4) 
# qrex <- q
# genestokeep <- subset( m.deseq, variable == 'issig' & value == 1 & meta == 'GFP.NEG',
#     select = 'GeneSymbol' , drop = TRUE)
# df2 <- subset( dftoplot, GeneSymbol %in% genestokeep) 
# q <- ggplot( df2, aes(GFP.NEG_log2Fold, HVneg.HVpos_2iLIF_avg_log2Fold)) 
# q <- f.scatterprettier(q, 140)
# q <- q + xlim(-2,4) 
# qhex <- q
# genestokeep <- subset( m.deseq, variable == 'issig' & value == 1 & meta == 'ALL.DIFF',
#     select = 'GeneSymbol' , drop = TRUE)
# df2 <- subset( dftoplot, GeneSymbol %in% genestokeep) 
# q <- ggplot( df2, aes(ALL.DIFF_log2Fold, `ESC.EpiLC:day1_log2Fold`)) 
# q <- f.scatterprettier(q,120)
# q <- q + xlim(-3,3)
# qepi <- q
# @

# <<>>=
# f.gridannotedscatter <- function(gplot, vp, posylab='hello',negylab='bye',
#     negxlab = '', posxlab = '', subtitle=''){
#    print(gplot, vp = vp)
#    downViewport(vp$name)
# #    longerstr <- ifelse(nchar(posylab)>=nchar(negylab), posylab, negylab)
#    vpdepth <- downViewport('axis-l.3-3-3-3')
#    xpos <- unit(-2.5, 'strwidth', 'A')
#    grid.text( posylab, x = xpos, y= 0.6, just=c('center','bottom'))
#    grid.text( negylab, x = xpos, y= 0.4, just=c('center','top'))
#     grid.lines( x=xpos,y=c(0.45,0.55),gp=gpar(fill='black',lwd=1),
# 	arrow=arrow(angle=45,type='closed', 
# 	    length=unit(0.5,'strwidth','A'), ends='both'))
#    upViewport(vpdepth)
#    vpdepth <- downViewport('axis-b.4-4-4-4')
#    ypos <- unit(-1.5, 'strheight', 'A')
#    grid.text( negxlab, x=0.3, y=ypos, just=c('right','center'))
#    grid.text( posxlab, x=0.7, y=ypos, just=c('left','center'))
#     grid.lines( x=c(0.4,0.6),y=ypos,gp=gpar(fill='black',lwd=1),
# 	arrow=arrow(angle=45,type='closed', 
# 	    length=unit(0.5,'strwidth','A'), ends='both'))
#    subtitleypos <- unit(-3.5, 'strheight','A')
#    grid.text( subtitle, x=0.5, y=subtitleypos, just=c('center','center'))
#    upViewport(vpdepth)
#    upViewport(1)
#    return(NULL)
# }
# f.drawrowscatter <- function(){
#     unit.oneplot <- unit.c( stringWidth('EpiLC+'),unit(1,'null'))
#     scattergrid <- grid.layout(3, 7, 
# 	heights = unit( c(1,1.5,1.5), c('null','strheight','strheight'),
# 	    list(NULL,'V','V')),
# 	widths = unit.c(unit(2,'strheight','A'), unit.oneplot, 
# 	    unit.oneplot, unit.oneplot 
# 	    )
#     )
#     pushViewport(viewport(layout=scattergrid))
#     vprex <- viewport(layout.pos.row=1,layout.pos.col=3,name='vprex')
#     vphex <- viewport(layout.pos.row=1,layout.pos.col=5,name='vphex')
#     vpepi <- viewport(layout.pos.row=1,layout.pos.col=7,name='vpepi')
#     vpylabel <- viewport(layout.pos.row=1,layout.pos.col=1,name='vpylabel')
#     f.gridannotedscatter( qrex, vprex, 'Rex1-', 'Rex1+','VNP+','VNP-')
#     f.gridannotedscatter( qhex, vphex, 'Hex+', 'Hex-','VNP+','VNP-',
# 	subtitle='log2( RNA fold change)')
#     f.gridannotedscatter( qepi, vpepi, 'EpiLC', 'ESC','Stem','Diff')
#     grid.text('log2( RNA fold change )', rot=90, vp=vpylabel)
#     upViewport(1)
# }
# @

# \section{Figure assembly}

# Assemble top:
# <<>>=
# f.labelsubplot <- function( subplotlabel, x=stringWidth('A'), y=1, ...) {
#     grid.text(subplotlabel,x=x,y=y, just=c('left','top'), 
# 	gp=gpar(fontsize=12, fontface='bold'), ...)
# }
# f.maketophalf <- function(){
#     pushViewport(viewport(layout = 
# 	grid.layout(1,2, widths = c(0.3,0.7))))
#     vpfold<-viewport(layout.pos.col=1,name='vpfold')
#     vpgo<-viewport(layout.pos.col=2,name='vpgo')
#     f.drawfoldchange(qfold,vp=vpfold)
#     pushViewport(vpgo)
#     f.drawGOanalysis()
#     upViewport(1)
#     f.labelsubplot('A', vp = vpfold)
#     f.labelsubplot('B', vp = vpgo)
#     upViewport(1)
# }
# f.makebottomhalf <- function(){
#     pushViewport(viewport(layout=grid.layout(1,3,widths=c(0.575,0.025,0.4))))
#     pushViewport( viewport(layout.pos.col=1))
#     f.drawrowscatter()
#     f.labelsubplot('C', y=0.95)
#     upViewport(1)
#     vptfs <- viewport(layout.pos.col=3)
#     print(qtfxdeseq, vp = vptfs)
#     f.labelsubplot('D', y=0.95, vp = vptfs)
#     upViewport(1)
# }
# f.plotfig2 <- function(){
#     pushViewport( viewport(gp=gpar(fontsize=base_size*0.8)))
#     pushViewport( viewport( layout=grid.layout(3,1,heights=c(6,0.25,6))))
#     pushViewport( viewport( layout.pos.row=1, name = 'vptop'))
#     f.maketophalf()
#     upViewport(1)
#     pushViewport( viewport( layout.pos.row=3, name = 'vptop'))
#     f.makebottomhalf()
#     upViewport(1)
# }
# @
# <<>>=
# setEPS()
# postscript('./outputdata/fig2/figHitBasics2.eps',width=7,height=4)
# f.plotfig2()
# dev.off()
# @


