dir.create(file.path('outputdata/figmacarthur'), showWarnings = FALSE)

library(reshape2)
library(plyr)
library(stringr)
library(ggplot2)
library(grid)

base_size <- 10



load('./intermediate_results/m.deseq.RData')
load('./intermediate_results/m.mac.RData')


dftoplot <- rbind(m.deseq, m.mac)
dftoplot <- subset(dftoplot, variable %in% c('log2Fold','issig'))
dftoplot <- dcast(dftoplot, 'GeneSymbol ~ meta + variable')
dftoplot[,-1] <- sapply(dftoplot[,-1], function(x) ifelse(is.finite(x),x,NA) )



f.scattermac <- function(gplotin, countlim = 60){
    gplotout <- gplotin + coord_fixed(ratio=2)
    gplotout <- gplotout + 
	geom_bin2d(binwidth=c(0.2,0.1)) 
    rain.colors <- rev( rainbow(7)[1:6])
    gplotout <- gplotout + scale_fill_gradientn(
	colours=rain.colors,
	breaks = c(0,countlim),
	guide = guide_colourbar(direction='horizontal',title=NULL,
	    label.position='top', barwidth = 2, barheight=0.25,
	    ticks = FALSE, )
    )
    gplotout <- gplotout + theme_grey(base_size=base_size)
    gplotout <- gplotout + theme(
	legend.position = c(0.80,0.12),
	legend.background = element_blank(),
# 	plot.margin=unit(rep(0,4),'lines'),
 	axis.text = element_text( colour='black'),
	axis.ticks = element_line(colour='black')
# 	axis.title.x = element_blank(),
# 	axis.title.y = element_blank()
# 	panel.grid.major = element_blank()
    )
    gplotout <- gplotout + scale_y_continuous(limits=c(-2,2))
#     gplotout <- gplotout + geom_hline(yintercept=0) + geom_vline(xintercept=0)
    gplotout <- gplotout 
    return(gplotout)
}


genestokeep <- subset( m.deseq, variable == 'issig' & value == 1 & 
    meta == 'GFP.NEG',  select = 'GeneSymbol' , drop = TRUE)
df2 <- subset( dftoplot, GeneSymbol %in% genestokeep) 
q <- ggplot( df2, aes(GFP.NEG_log2Fold, t36h_log2Fold))  
q <- f.scattermac(q,70) + xlim(-2.5,5)
q36 <- q
q <- ggplot( df2, aes(GFP.NEG_log2Fold, t84h_log2Fold))  
q <- f.scattermac(q,40) + xlim(-2.5,5)
q84 <- q
q <- ggplot( df2, aes(GFP.NEG_log2Fold, t132h_log2Fold))  
q <- f.scattermac(q,40) + xlim(-2.5,5)
q132 <- q



genestokeep <- subset( m.deseq, variable == 'issig' & value == 1 & 
    meta == 'ALL.DIFF',  select = 'GeneSymbol' , drop = TRUE)
df2 <- subset( dftoplot, GeneSymbol %in% genestokeep) 
q <- ggplot( df2, aes(ALL.DIFF_log2Fold, t36h_log2Fold))  
q <- f.scattermac(q,90) + xlim(-3.75, 3.75)
q36a <- q
q <- ggplot( df2, aes(ALL.DIFF_log2Fold, t84h_log2Fold))  
q <- f.scattermac(q,70) + xlim(-3.75, 3.75)
q84a <- q
q <- ggplot( df2, aes(ALL.DIFF_log2Fold, t132h_log2Fold))  
q <- f.scattermac(q,50) + xlim(-3.75, 3.75)
q132a <- q



setEPS()
postscript('./outputdata/figmacarthur/scatter.mac.eps',width=7.5,height=5)
pushViewport( viewport( layout = grid.layout( 2, 3 ) ))
print(q36, vp = viewport(layout.pos.col=1,layout.pos.row=1))
print(q84, vp = viewport(layout.pos.col=2,layout.pos.row=1))
print(q132, vp = viewport(layout.pos.col=3,layout.pos.row=1))
print(q36a, vp = viewport(layout.pos.col=1,layout.pos.row=2))
print(q84a, vp = viewport(layout.pos.col=2,layout.pos.row=2))
print(q132a, vp = viewport(layout.pos.col=3,layout.pos.row=2))
dev.off()


