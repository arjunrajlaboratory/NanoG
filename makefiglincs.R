dir.create(file.path('outputdata/figlincs'), showWarnings = FALSE)

library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(grid)

Stangle('s_functions.Rnw');source('s_functions.R')
Stangle('funcs_figsource.Rnw');source('funcs_figsource.R')


load('./intermediate_results/m.lincs.RData')
load('./intermediate_results/updatelincs.RData')

base_size <- 10


# That shows that the new annotations are consistent with the old ones from the
# paper. However, many (71) old genes are now missing. 
# Of the 166 new listed genes, 151 have a correspondence with the old.  

# The following function transforms a list of old names into new names.

f.newlinc.fromoldlinc <- function( oldnames.in, updatelincs){
    #updatelincs has columns oldname and newname
    newlinc <- with(updatelincs, newname[ oldname %in% oldnames.in] )
    newlinc <- unique(as.character( newlinc))
    newlinc <- newlinc[ !is.na(newlinc) ]
    return(newlinc)
}


# section{Basic stats on the lincs}
# Basic stats on the lincs:


m.lincsonly <- subset(m.lincs, grepl('^linc[[:digit:]]', GeneSymbol))



dftoplot <- subset(m.lincsonly, meta %in% c('ALL.DIFF','GFP.NEG'))
dftoplot <- dcast( dftoplot, 'GeneSymbol ~ meta + variable' )
print( table( dftoplot[,c('GFP.NEG_signedsig','ALL.DIFF_signedsig')]))
print( addmargins( table( dftoplot[, c('GFP.NEG_issig','ALL.DIFF_issig')])))



dftoplot <- subset(m.lincs, grepl('^linc[[:digit:]]', GeneSymbol))
dftoplot <- subset(dftoplot, meta %in% c('ALL.DIFF','GFP.NEG'))
dftoplot <- dcast(dftoplot, 'GeneSymbol ~ meta + variable')
q <- ggplot( subset(dftoplot, GFP.NEG_issig & ALL.DIFF_issig), 
    aes( x = GFP.NEG_log2Fold, y = ALL.DIFF_log2Fold ) )
q <- q + geom_point()
q <- q + xlim(-2,4) + ylim(-3,3)
qscatter <- q + ggtitle('sig in both')
qeither <- q %+% subset(dftoplot, GFP.NEG_issig | ALL.DIFF_issig)
qeither <- qeither + ggtitle('sig in either')
qall <- q %+% dftoplot
qall <- qall + ggtitle('all ES lincs')



dftoplot <- transform(dftoplot, GFP.NEG_issig = factor(GFP.NEG_issig, c(1,0)))
q <- ggplot( dftoplot, aes( x = GFP.NEG_log2Fold, fill = GFP.NEG_issig ) )
q <- q + geom_bar(binwidth=0.33)
q <- q + scale_fill_manual(values=c('black','grey'), guide=FALSE) 
q <- q + theme_minimal() + xlim(-2,4)
qgfpneg <- q
dftoplot <- transform(dftoplot, ALL.DIFF_issig = factor(ALL.DIFF_issig, c(1,0)))
q <- ggplot( dftoplot, aes( x = ALL.DIFF_log2Fold, fill = ALL.DIFF_issig ) )
q <- q + geom_bar(binwidth=0.33)
q <- q + scale_fill_manual(values=c('black','grey'), guide=FALSE) 
q <- q + theme_minimal() + xlim(-3,3)
qalldiff <- q



pdf('./outputdata/figlincs/lincoverview.pdf',width=9,height=3)
pushViewport(viewport(layout=grid.layout(1,3)))
print( qgfpneg, vp = viewport(layout.pos.col=1))
print( qalldiff, vp = viewport(layout.pos.col=2))
print( qscatter, vp = viewport(layout.pos.col=3))
dev.off()



pdf('./outputdata/figlincs/lincscatters.pdf',width=9,height=3)
pushViewport(viewport(layout=grid.layout(1,3)))
print( qall, vp = viewport(layout.pos.col=1))
print( qeither, vp = viewport(layout.pos.col=2))
print( qscatter, vp = viewport(layout.pos.col=3))
dev.off()


# section{Pluripotency-implicated lincs}
# Load the Nanog luciferase assay data from Guttman 2011

luciflinc <- read.csv('./externaldata/guttman11/supp7.csv')
luciflinc <- subset(luciflinc, Nanog_Z1 < -6 & Nanog_Z2 < -6 
    & grepl('linc',Knockdown) )
luciflinc <- melt(luciflinc, measure.vars = c('Alamar_Puro','Alamar_Puro2'))
q <- ggplot( luciflinc, aes(value/Alamar_NoPuro, Knockdown, color = variable))
q <- q + geom_point()
print(q)


# The plot shows that linc1463 is the outlier in terms of low viability. Discard it.

luciflinc <- subset( luciflinc, value/Alamar_NoPuro > 0.6 )
plurilincs <- unique(as.character(luciflinc$Knockdown))


# The result is 26 lincs, wich is the same number they quote in the paper.
# Let's look at these genes


dftoplot <- subset( m.lincs, meta %in% c('GFP.NEG', 'ALL.DIFF'))
dftoplot <- dcast( dftoplot, 'GeneSymbol + meta ~ variable' )
dftoplot <- transform(dftoplot, meta = factor( meta, c('GFP.NEG','ALL.DIFF')))



plurilincsnew <- f.newlinc.fromoldlinc( plurilincs, updatelincs )



q <- f.logfold.dotplot( subset(dftoplot, GeneSymbol %in% plurilincsnew))
q <- q + xlim(-2,3) + ggtitle('Pluripotency-associated')
setEPS()
postscript('./outputdata/figlincs/plurilincs.dotplot.eps',width=2.4,height=4)
print(q)
dev.off()


# section{lincs that repress lineages}

# These are read directly off of Figure 3

linlincs <- list()
linlincs[['Ectoderm']] <- c('1253','1331','1557','1562','1602','1230')
linlincs[['Neuroectoderm']] <- c('1230','1335','1463','1582','1612')
linlincs[['Endoderm']] <- c('1242','1304','1356','1400','1456','1490',
    '1517','1588','1617','1623','1627','1633','1388','1390')
linlincs[['Mesoderm']] <- c('1388','1390','1307','1434','1526','1552')
linlincs[['Trophectoderm']] <- c('1470','1604')
linlincs <- melt(linlincs)
names(linlincs) <- c('GeneSymbol','lineage')
linlincs$GeneSymbol <- paste0('linc',linlincs$GeneSymbol)
linlincs <- ddply(linlincs, 'GeneSymbol', summarize, 
    lineage = paste( lineage, collapse = ','))
# shorten names when a linc is assigned to two lineages
linlincs <- transform( linlincs, 
    lineage = str_replace( lineage, '^(....).*,\\s*(....).*', '\\1, \\2'))
linorder <- c('Ectoderm', 'Ecto, Neur', 'Neuroectoderm','Mesoderm',
    'Endo, Meso', 'Endoderm', 'Trophectoderm')
linlincs <- transform( linlincs, lineage = factor(lineage, linorder))
linlincsnew <- ddply( linlincs, 'lineage', summarize,
    GeneSymbol = f.newlinc.fromoldlinc(GeneSymbol, updatelincs))



df2 <- merge(dftoplot, linlincsnew, by='GeneSymbol')
q <- f.logfold.dotplot(df2) + xlim(-2,3)
q <- q + facet_grid('lineage ~.',scales='free_y',space='free_y')
q <- q + theme( strip.text.y = element_text(angle = 0))
qrepressors <- q
q <- q + ggtitle('linc repressors of lineages')
setEPS()
postscript('./outputdata/figlincs/linrepressors.dotplot.eps',width=3.3,height=4)
print(q)
dev.off()


# section{Make a figure}


setEPS()
postscript('./outputdata/figlincs/figlinc.eps',width=5,height=3)
pushViewport(viewport(layout=grid.layout(2,2,widths=c(1,1.8))))
print( qgfpneg, vp = viewport(layout.pos.col=1,layout.pos.row=1))
print( qalldiff, vp = viewport(layout.pos.col=1,layout.pos.row=2))
print( qrepressors, vp = viewport(layout.pos.col=2, layout.pos.row=c(1,2)))
dev.off()




