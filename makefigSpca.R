dir.create(file.path('outputdata/figSpca'), showWarnings = FALSE)

# \section{PCA based overview of the data}

# \subsection{Variance Stablization and PCA computation}

# Following DESeq vignette

library(DESeq)
library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(grid)
Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')


base_size <- 10

load('./intermediate_results/cds.RData')
cdsBlind <- DESeq::estimateDispersions( cds, method='blind')



vsd <- DESeq::getVarianceStabilizedData( cdsBlind )
# Remove difference between A and C samples here: 
# vectortoremove <- as.matrix(c(1,1,1,1,1,-1,-1,-1,-1,-1))/sqrt(10)
# row.names(vectortoremove) <- c('A1','A2','A3','A4','A5',
# 	'C1','C2','C3','C4','C5')
# vsd <- vsd - vsd %*% vectortoremove %*% t(vectortoremove)

# Remove the average of all samples
vectortoremove <- as.matrix(c(1,1,1,1,1,1,1,1,1,1))/sqrt(10)
row.names(vectortoremove) <- c('A1','A2','A3','A4','A5',
	'C1','C2','C3','C4','C5')
vsd <- vsd - vsd %*% vectortoremove %*% t(vectortoremove)
# Change column names and melt.
m.vsd <- f.meltrumcounts( vsd ) 
m.vsd$variable <- 'DESeqVST'
m.vsd$source <- 'deseq'
m.vsd <- m.vsd[,c('source','meta','GeneSymbol','variable','value')]

vsdmat <- dcast(m.vsd, 'GeneSymbol ~ meta')
vsdmat <- as.matrix(vsdmat[,-1])
pca <- prcomp( vsdmat, center=TRUE)
m.pca <- as.data.frame(pca$rotation)
# print(pca$rotation[,c(1:4)])
# print(summary(pca))
m.pca <- cbind( meta = row.names(m.pca), m.pca)
m.pca <- f.colsplitandmerge(m.pca,'meta','.',c('cond','rep'))
m.pca$cond <- factor(f.sample.labeller( value = m.pca$cond ))


q <- ggplot( m.pca, aes( x=PC1, y=PC2, color=cond, shape=rep))
q <- q + geom_point() + coord_equal()
q <- q + theme_minimal(base_size =base_size)
q <- q + guides(
    color = guide_legend(title='sample', ncol=1, byrow=TRUE,
	keyheight = stringHeight('D')),
    shape = guide_legend(title='replicate',direction='horizontal',
	title.position='top',keyheight=stringHeight('D'))
)
qtoppca <- q


pcavars <- data.frame( pcname = colnames(pca$rotation),variance = pca$sdev^2 )
pcavars$pcname <- with(pcavars, factor(pcname,pcname))
q <- ggplot( pcavars, aes(pcname, variance)) + geom_bar(stat='identity')
q <- q + theme_minimal(base_size=base_size)
q <- q + theme( 
    axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
q <- q + xlab(NULL)
qpcavar <- q


dftoplot <- m.pca
q <- ggplot(m.pca, aes(y= factor(cond, levels=rev(levels(cond))), x=PC3, 
    shape=rep,color=cond)) 
q <- q + geom_vline(xintercept=0) + geom_point()
q <- q + theme_minimal(base_size=base_size) 
q <- q + theme(
    legend.position='none',
    plot.margin = unit(rep(0,4),'lines'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.margin = unit(0,'lines'),
    axis.text.y = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
q <- q + ylab(NULL)
qpc3 <- q 



f.drawvars.andpc3 <- function(){
print(qpcavar, vp = current.viewport())
vpdepth <- downViewport('panel.3-4-3-4')
vpinsert <- viewport(x=0.65,y=0.65,width=0.7,height=0.7,
    clip='off',name='inset',gp=gpar(fontsize=base_size))
grid.rect( gp = gpar(lty=0,fill='white'), vp=vpinsert)
print(qpc3, vp=vpinsert)
upViewport(vpdepth)
}



f.labelsubplot <- function( subplotlabel, x=stringWidth('A'), y=1, ...) {
    grid.text(subplotlabel,x=x,y=y, just=c('left','top'), 
	gp=gpar(fontsize=12, fontface='bold'), ...)
}
f.drawpcafig <- function(){
pushViewport( viewport( layout = grid.layout(1,2,widths=c(5,3))))
vppcavar <- viewport( layout.pos.col = 2, name='vppcavar')
vptoppca <- viewport( layout.pos.col = 1, name = 'vptoppca')
print(qtoppca, vp = vptoppca)
f.labelsubplot('A', vp = vptoppca)
pushViewport(vppcavar)
f.drawvars.andpc3()
upViewport(1)
f.labelsubplot('B', vp = vppcavar)
}


setEPS()
postscript('./outputdata/figSpca/figSpca.eps',width=6,height=2)
f.drawpcafig()
dev.off()

