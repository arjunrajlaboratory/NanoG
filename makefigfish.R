dir.create(file.path('outputdata/figfish'), showWarnings = FALSE)

library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)
pathToFISH <- 'inputdata/'


# \section{Nanog, Crabp2 and T}

outputTextFile <- './outputdata/figfish/countsForVennDiagrams.txt'
capture.output(print('counts for Venn diagrams'), file = outputTextFile, append = FALSE)

printToStatsFile <- function(x){
    capture.output(print(x),
        file = outputTextFile,
        append = TRUE)
}

datnew <- list()
datnew[[1]] <- cbind( read.csv(paste0(pathToFISH,'new130727.csv')), date = '130727')
datnew[[2]] <- cbind( read.csv(paste0(pathToFISH,'new130729.csv')), date = '130729')
datnew <- do.call(rbind,datnew)
datnew <- subset(datnew, isGood == 1)



Nanogthresh <- 30
Tthresh <- 50
Crabp2thresh <- 30
datnew <- transform(datnew, 
 	Nanog.Lev = factor(ifelse(Nanog.RNA>Nanogthresh, 'NanogHigh','NanogLow')),
 	T.Lev = factor(ifelse(T.RNA>Tthresh, 'THigh', 'TLow')), 
 	Crabp2.Lev = factor(ifelse(Crabp2.RNA>Crabp2thresh, 'Crabp2High', 'Crabp2Low')) )





dftoplot <- melt(datnew, measure.vars =c('Nanog.RNA','T.RNA','Crabp2.RNA'))
dftoplot <- subset( dftoplot, variable %in% c('T.RNA','Crabp2.RNA'))
q <- ggplot(dftoplot, aes(value)) + geom_dotplot()
q <- q + facet_grid('Nanog.Lev~variable',scales='free_x')
q <- q + theme_minimal()
setEPS()
postscript('./outputdata/figfish/T.Crabp2.eps',width=3.8,height=3)
print(q)
dev.off()


printToStatsFile(table(dftoplot[,c('T.Lev','Crabp2.Lev','Nanog.Lev')]))



q <- ggplot(datnew, aes(Nanog.RNA)) + geom_histogram(binwidth=10)
q <- q + geom_vline(x = Nanogthresh) + theme_minimal()
setEPS()
postscript('./outputdata/figfish/NanogDist.eps', width =3.8, height = 1.9)
print(q)
dev.off()


# section{Nanog, Rex, VNP}


datrexvnpnew <- read.csv(paste0(pathToFISH,'new130821.csv'))
datrexvnpnew <- subset(datrexvnpnew, isGood==1)
VNPthresh <- 40
Nanogthresh <- 30
Rexthresh <- 50
datrexvnpnew <- transform(datrexvnpnew, 
	VNP.Lev = ifelse(VNP.RNA<VNPthresh, 
		paste0('VNP','<',VNPthresh),
		paste0('VNP','>',VNPthresh)),
	Nanog.Lev = ifelse(Nanog.RNA<Nanogthresh, 
		paste0('Nanog','<',Nanogthresh),
		paste0('Nanog','>',Nanogthresh)),
	Rex1.Lev = ifelse(Rex1.RNA<Rexthresh, 
		paste0('Rex1','<',Rexthresh),
		paste0('Rex1','>',Rexthresh))
)




base_size <- 10
f.fishjointdensity <- function(gplotin, countlim = 60){
    gplotout <- gplotin + geom_bin2d(binwidth=c(30,30),color='black')
    gplotout <- gplotout + coord_equal()
    gplotout <- gplotout + theme_grey(base_size=base_size)
    rain.colors <- rev( rainbow(7)[1:6])
    gplotout <- gplotout + scale_fill_gradientn(
	colours=rain.colors,
	limits = c(1,95),
	breaks = c(0, countlim),
	trans='sqrt',
	guide = guide_colourbar(
	    direction='vertical',title='#cells',
	    label.position='right', barwidth =0.5, 
	    barheight = 5, ticks = FALSE
	)
    )
    return(gplotout)
}
q <- ggplot(datrexvnpnew, aes(Nanog.RNA,VNP.RNA))  
q <- f.fishjointdensity(q,countlim=c(10,30,50,80))
q <- q + xlab('Nanog RNA') + ylab('VNP RNA')
qvnp <- q + ylim(0, 675)
q <- ggplot(datrexvnpnew, aes(Nanog.RNA,Rex1.RNA))  
q <- f.fishjointdensity(q,countlim=c(10,30,50,80))
q <- q + xlab('Nanog RNA') + ylab('Rex1 RNA')
qrex <- q + ylim(0, 675)

f.labelsubplot <- function( subplotlabel, x=stringWidth('A'), y=1, ...) {
    grid.text(subplotlabel,x=x,y=y, just=c('left','top'), 
	gp=gpar(fontsize=12, fontface='bold'), ...)
}
f.drawfishcorr <- function(){
    pushViewport(viewport(layout=grid.layout(1,2)))
    vpvnp <- viewport(layout.pos.col=1)
    vprex <- viewport(layout.pos.col=2)
    print(qvnp, vp = vpvnp)
    print(qrex, vp = vprex)
    f.labelsubplot('A',vp=vpvnp)
    f.labelsubplot('B',vp=vprex)
    upViewport()
}
require(grid)
setEPS()
postscript('./outputdata/figfish/rex.vnp.Dist.eps',width=6,height=2.8)
f.drawfishcorr()
dev.off()




printToStatsFile(addmargins(table(datrexvnpnew[,c('VNP.Lev','Nanog.Lev')])))

nanogAndVNPMeans <- ddply(datrexvnpnew, c('VNP.Lev','Nanog.Lev'), summarize, 
	meanNanog.RNA = mean(Nanog.RNA), meanVNP.RNA = mean(VNP.RNA))

printToStatsFile(nanogAndVNPMeans)

printToStatsFile(addmargins(table(datrexvnpnew[,c('Rex1.Lev','Nanog.Lev')])))


printToStatsFile('for Rex, VNP fold change:')
printToStatsFile(ddply(datrexvnpnew, 'Nanog.Lev', summarize, Rex1.RNA=mean(Rex1.RNA),
    VNP.RNA=mean(VNP.RNA)))

# section{Oct4, T, Tbx6}


datnew <- read.csv(paste0(pathToFISH,'new130903.csv'))
datnew <- subset(datnew, isGood == 1)
datnew <- transform(datnew, datfilenum = as.numeric(str_match(dataFile,'[[:digit:]]+')))
# Shows that there is a need to exclude Oct4 RNA values from the later data files.
# ggplot(datnew, aes(datfilenum, Oct4.RNA)) + geom_point()
# datnew <- transform(datnew, Oct4.RNA = ifelse(datfilenum > 74, NA, Oct4.RNA))
datnew <- subset(datnew, datfilenum <= 74)



datnew <- transform(datnew, Oct4.Lev = factor(Oct4.RNA > 80, c(TRUE,FALSE),
    labels = c('Oct4>80', 'Oct4<80')))
datnew <- transform(datnew, T.Lev = factor(T.RNA > 50, c(FALSE,TRUE),
    labels = c('T<50', 'T>50')))
datnew <- transform(datnew, Tbx6.Lev = factor(Tbx6.RNA > 10, c(TRUE,FALSE),
    labels = c('Tbx6>10', 'Tbx6<10')))



printToStatsFile(table(datnew[,c('T.Lev','Tbx6.Lev','Oct4.Lev')]))


q <- ggplot(datnew, aes(T.RNA, Tbx6.RNA, color = Oct4.Lev)) + geom_jitter()
q <- q + guides(color=guide_legend(title=NULL,direction='horizontal' ))
q <- q + theme_minimal()
q <- q + theme(legend.position = 'bottom')
q <- q + ylim(-1, 30)
setEPS()
postscript('./outputdata/figfish/T.Tbx6.byOct4.Dist.eps', width=3,height=3)
print(q)
dev.off()



dftoplot <- melt(datnew, measure.vars =c('Oct4.RNA','T.RNA','Tbx6.RNA'))
dftoplot <- subset( dftoplot, variable %in% c('T.RNA','Tbx6.RNA'))
q <- ggplot(dftoplot, aes(value)) + geom_dotplot()
q <- q + facet_grid('Oct4.Lev~variable',scales='free_x')
q <- q + theme_minimal()
setEPS()
postscript('./outputdata/figfish/T.Tbx6.eps',width=3.8,height=3)
print(q)
dev.off()



q <- ggplot(datnew, aes(Oct4.RNA)) + geom_histogram(binwidth=20)
q <- q + geom_vline(x = 80) + theme_minimal()
setEPS()
postscript('./outputdata/figfish/Oct4Dist.eps', width =3.8, height = 1.9)
print(q)
dev.off()


# section{Oct4, Nanog}


datnew <- read.csv(paste0(pathToFISH,'new131008.csv'))
datnew <- subset(datnew, isGood == 1)
datnew <- transform(datnew, Nanog.Lev = factor(Nanog.RNA>30, c(TRUE,FALSE),
    labels = c('Nanog > 30', 'Nanog < 30')))
datnew <- transform(datnew, Oct4.Lev = factor(Oct4.RNA > 80, c(TRUE,FALSE),
    labels = c('Oct4 > 80', 'Oct4 < 80')))


printToStatsFile(table(datnew[,c('Oct4.Lev','Nanog.Lev')]))

printToStatsFile(addmargins(table(datnew$Nanog.Lev)))

printToStatsFile(addmargins(prop.table(table(datnew$Nanog.Lev))))

printToStatsFile('for Oct4 fold change:')
printToStatsFile(ddply(datnew, 'Nanog.Lev', summarize, Oct4=mean(Oct4.RNA)))


q <- ggplot(datnew, aes(Nanog.RNA)) + geom_histogram(binwidth=10)
q <- q + geom_vline(x = 30) + theme_minimal()
qNanog <- q
q <- ggplot(datnew, aes(Oct4.RNA)) + geom_histogram(binwidth=20)
q <- q + geom_vline(x = 80) + theme_minimal()
qOct4 <- q
f.drawNanogOct4 <- function(){
    pushViewport(viewport(layout=grid.layout(1,2)))
    vpNanog <- viewport(layout.pos.col=1)
    vpOct4 <- viewport(layout.pos.col=2)
    print(qNanog, vp = vpNanog)
    print(qOct4, vp = vpOct4)
    f.labelsubplot('A',vp=vpNanog)
    f.labelsubplot('B',vp=vpOct4)
    upViewport()
}
setEPS()
postscript('./outputdata/figfish/NanogOct4Dist.eps', width =6, height = 3)
f.drawNanogOct4()
dev.off()

