### R code from vignette source 's_functions.Rnw'

###################################################
### code chunk number 1: s_functions.Rnw:4-15
###################################################
f.plot_counts_alacarte<-function(vectorofgenes, countlist, condlist, ...) {
    for (i in 1:length(vectorofgenes)) {
	datatoplot<-countlist[ vectorofgenes[i], ]
	max_y<-max(datatoplot)
	plot.default( condlist,
             datatoplot, xlab = 'Condition', ylab = 'counts',
             ylim = c(0,max_y), main = vectorofgenes[i], xaxt = 'n', ...)
        axis(1, at = seq(along = levels(condlist) ),
             labels = levels( condlist ), las = 2 )
    }
}


###################################################
### code chunk number 2: s_functions.Rnw:22-31
###################################################
f.deseq_do_andwrite <- function(countdataset, condpair,
                                fileprefix ) {
    message( paste( 'Testing', condpair[1], 'vs', condpair[2]))
    nbinomresult <- nbinomTest(countdataset, condpair[1], condpair[2])
    filename <- paste0(fileprefix, condpair[1] , condpair[2])
    write.table( nbinomresult, file = filename, quote = FALSE, sep = '\t')
    return( list(condA = condpair[1], condB = condpair[2],
                 filename = filename))
}


###################################################
### code chunk number 3: s_functions.Rnw:35-41
###################################################
f.deseq_testallpairs <- function(countdataset, fileprefix ) {
    FUN = function (x)  f.deseq_do_andwrite(countdataset,x,fileprefix)
    nbinomindex<- apply(X = combn(levels(conditions(countdataset)),m=2),
                        MARGIN = 2, FUN = FUN)
    return(nbinomindex)
}


###################################################
### code chunk number 4: s_functions.Rnw:50-68
###################################################
f.loadntest<-function(condA,condB,nindex){
    ntable<-t(sapply(nindex,
                     FUN=function (x) c(paste0(x$condA,'.',x$condB),
                     x$filename)))
    forward.name<-paste0(condA,'.',condB)
    reverse.name<-paste0(condB,'.',condA)
    if (forward.name %in% ntable[,1] ) {
        dframe<-read.delim(ntable[match(forward.name,ntable[,1]),2])}
    else if (reverse.name %in% ntable[,1]) {
        dframe<-read.delim(ntable[match(reverse.name,ntable[,1]),2])
        dframe[,c('baseMeanA','baseMeanB')]<-
            dframe[,c('baseMeanB','baseMeanA')]
        dframe$foldChange<-1/dframe$foldChange
        dframe$log2FoldChange<- - dframe$log2FoldChange
        }
    else dframe=NULL
    return(dframe)
}


###################################################
### code chunk number 5: s_functions.Rnw:76-97
###################################################
nbinomTestSimilarity <- function(cds, condA, condB, maxratiovec = c(2) ){
    colA <- conditions(cds) == condA
    colB <- conditions(cds) == condB
    rawScvA <- fData(cds)[,paste('disp',dispTable(cds)[condA],sep='_')]
    rawScvB <- fData(cds)[,paste('disp',dispTable(cds)[condB],sep='_')]
    outdf <- data.frame( id = rownames(counts(cds)), stringsAsFactors=FALSE)
    for ( maxratio in maxratiovec) {
	pval1 <- nbinomTestForMatricesSimilarity( counts(cds)[,colA],
	    counts(cds)[,colB], sizeFactors(cds)[colA], sizeFactors(cds)[colB],
	    rawScvA, rawScvB, maxratio)
	pval2 <- nbinomTestForMatricesSimilarity( counts(cds)[,colB],
	    counts(cds)[,colA], sizeFactors(cds)[colB], sizeFactors(cds)[colA],
	    rawScvB, rawScvA, maxratio)
	# combine the p values as though the two tests are independent:
	pval <- 1 - (1- pval1) * (1-pval2)
	padj <- p.adjust( pval, method = 'BH')
	outdf[[ paste0('simpval_', as.character(maxratio)) ]] <- pval
	outdf[[ paste0('simpadj_', as.character(maxratio)) ]] <- padj
    }
    return( outdf )
}


###################################################
### code chunk number 6: s_functions.Rnw:100-135
###################################################
nbinomTestForMatricesSimilarity <- function(countsA, countsB, 
    sizeFactorsA, sizeFactorsB, dispsA, dispsB, ABratio){
    
    kAs <- rowSums( cbind( countsA ))
    kBs <- rowSums( cbind( countsB ))
    mus <- rowMeans( cbind( t(t(countsA)/sizeFactorsA), t(t(countsB)/sizeFactorsB)))
    muBs <- mus*(sum(sizeFactorsA)+sum(sizeFactorsB)) /
	(ABratio * sum(sizeFactorsA) + sum(sizeFactorsB) )
    muAs <- ABratio * muBs
    fullVarsA <- pmax(
	muAs * sum(sizeFactorsA) + dispsA * muAs^2 * sum(sizeFactorsA^2),
	muAs * sum(sizeFactorsA) * (1 + 1e-08))
    fullVarsB <- pmax(
	muBs * sum(sizeFactorsB) + dispsB * muBs^2 * sum(sizeFactorsB^2),
	muBs * sum(sizeFactorsB) * (1 + 1e-08))
    sumDispsA <- (fullVarsA - muAs * sum(sizeFactorsA))/(muAs*sum(sizeFactorsA))^2
    sumDispsB <- (fullVarsB - muBs * sum(sizeFactorsB))/(muBs*sum(sizeFactorsB))^2
    
    testonegene <- function (i) {
	if (kAs[i] == 0 & kBs[i] == 0) 
	   return(NA) 
	ks <- 0:(kAs[i] + kBs[i])
	ps <- dnbinom(ks, mu = muAs[i]*sum(sizeFactorsA), size = 1/sumDispsA[i]) *
	    dnbinom( kAs[i] + kBs[i] - ks, mu = muBs[i]*sum(sizeFactorsB),
		size = 1/sumDispsB[i])
	pobs <- dnbinom(kAs[i], mu = muAs[i]*sum(sizeFactorsA), size = 1/sumDispsA[i]) *
	    dnbinom( kBs[i], mu = muBs[i]*sum(sizeFactorsB),
		size = 1/sumDispsB[i])
	stopifnot(pobs==ps[kAs[i] + 1])
	numer <- ps[1:(kAs[i] + 1)]
	return( sum(numer)/sum(ps) )
    }

    return( sapply(seq(along=kAs), testonegene))
}


###################################################
### code chunk number 7: s_functions.Rnw:138-151
###################################################
f.loadsimtest <- function(condA,condB,nindex){
    ntable<-t(sapply(nindex,
                     FUN=function (x) c(paste0(x$condA,x$condB),
                     x$filename)))
    forward.name <- paste0(condA,condB)
    reverse.name <- paste0(condB,condA)
    if (forward.name %in% ntable[,1] ) 
        dframe <- read.delim(ntable[match(forward.name,ntable[,1]),2])
    else if (reverse.name %in% ntable[,1]) 
        dframe <- read.delim(ntable[match(reverse.name,ntable[,1]),2])
    else dframe = NULL
    return(dframe)
}


###################################################
### code chunk number 8: s_functions.Rnw:159-170
###################################################
f.summarize1v1 <- function( nbindex, countdataset, FUN = NULL){
    lconds <- levels(conditions(cds))
    summarymatrix <-
        array( data = 0, dim = replicate(2, length(lconds)),
          dimnames = replicate(2, lconds, simplify = FALSE))
    for (nbtest in nbindex) {
        dframe<-read.delim( nbtest$filename, header=TRUE)
        summarymatrix[ nbtest$condA, nbtest$condB] <- FUN(dframe)
     }
    return(summarymatrix)
}


###################################################
### code chunk number 9: s_functions.Rnw:176-181
###################################################
f.meetsfdr <- function (x, fdr = 0.1) {
    meetsfdr<- x$padj < fdr
    meetsfdr[is.na(meetsfdr)]<-FALSE
    return(meetsfdr)
}


###################################################
### code chunk number 10: s_functions.Rnw:185-188
###################################################
f.hitsatfdr <- function (x, fdr = 0.1) {
    sum( f.meetsfdr(x,fdr=fdr))
}


###################################################
### code chunk number 11: s_functions.Rnw:192-194
###################################################
f.isexpressed <- function (x) x$baseMean > 0
f.getlogbaseMean <- function (x) log10(x$baseMean)


###################################################
### code chunk number 12: s_functions.Rnw:205-222
###################################################
f.plothitstats <- function(condA, condB, nindex, fdr=0.1, regularizelog=0.1,
                           logmeanbreaks="FD", log2foldbreaks="FD") {
    # Plots log2foldchange and baseMean hit and all ditributions
    dframe <- f.loadntest(condA,condB,nindex)
    hitlist <- which(dframe$padj<fdr)
    meanhist<-hist( log10(regularizelog + dframe$baseMean),
                   breaks=logmeanbreaks, main=paste0(condA,'vs',condB),
                   xlab='log10 baseMean express.', border=NA,col='black')
    hist( log10(regularizelog + dframe$baseMean[hitlist]),
         breaks=meanhist$breaks, col='red',add=TRUE,border=NA)
    foldhist<- hist( dframe$log2FoldChange,
                   breaks=log2foldbreaks, main=paste0(condA,'vs',condB),
                   xlab='log2FoldChange', border=NA, col='black',
                    xlim=c(-3,3))
    hist(dframe$log2FoldChange[hitlist], breaks=foldhist$breaks,
         col='red',add=TRUE, border=NA)
}


###################################################
### code chunk number 13: s_functions.Rnw:227-245
###################################################
f.table_mostextreme <- function (dframe, nshow=10, ...) {
    biggest <- dframe[order(dframe$log2FoldChange,decreasing=TRUE)[1:nshow] ,
                      c("id","log2FoldChange","baseMean")]
    smallest <- dframe[order(dframe$log2FoldChange,decreasing=FALSE)[1:nshow] ,
                      c("id","log2FoldChange","baseMean")]
    names(biggest)<-c('Gene','log2Fold','MeanExpr')
    names(smallest)<-names(biggest)
    return( xtable( cbind(biggest, smallest), align="llcc|lcc", ... ))
    }
f.table_filteredextremes <- function(condA, condB, nindex, nshow=10,
                                     filterfunc= function (x) TRUE,
                                     label, caption, ...) {
    # filterfunc returns a logical used to index the rows of the dataframe
    dframe <- f.loadntest(condA,condB,nindex)
    dframe <- dframe[filterfunc(dframe),]
    xt <- f.table_mostextreme(dframe, nshow, caption=caption, label=label, ...)
    print(xt, caption.placement='top',include.rownames=FALSE, ...)
    }


###################################################
### code chunk number 14: s_functions.Rnw:253-267
###################################################
f.vectorsforgoseq <- function(condA, condB, nindex,
                              universefunc = function(x) rep(TRUE,dim(x)[1]),
                              isDEfunc = function(x) x$padj < 0.1,
                              biasfunc = function(x) NULL, ...) {
    # universefunc and isDEfunc should return logicals for every entry
    # in the DEseq results. biasfunc should return null or a numeric vector.
    dframe <- f.loadntest(condA,condB,nindex)
    isinuniverse <- universefunc(dframe)
    isDE <- isDEfunc(dframe)
    DEgenes <- isDE[isinuniverse]
    names(DEgenes) <- dframe$id[isinuniverse]
    bias.data <- biasfunc(dframe, ...)[isinuniverse]
    return( list(DEgenes=DEgenes, bias.data=bias.data) )
}


###################################################
### code chunk number 15: s_functions.Rnw:271-281
###################################################
f.goseqresults <- function( goseqvectors , genome='mm9', id='geneSymbol', ...) {
    pwf <- nullp(DEgenes = goseqvectors$DEgenes, genome=genome,
                 id=id, bias.data = goseqvectors$bias.data)
    GOdat<-goseq(pwf,genome=genome, id=id, ...)
    #GOdat$GOTERM<-sapply(mget(GOdat$category,GOTERM),function (x) x@Term)
    GOdat$GOTERM <- Term(GOdat$category)
    GOdat$over_padj<-p.adjust(GOdat$over_represented_pvalue,method="BH")
    GOdat$under_padj<-p.adjust(GOdat$under_represented_pvalue,method="BH")
    return(GOdat)
}


###################################################
### code chunk number 16: s_functions.Rnw:288-304
###################################################
f.getrpkmstats <- function(cds,cond,length_bases){
    # Assumes sizeFactors(cds) are millions of total reads
    col <- conditions(cds)==cond
    disps <- fData(cds)[, paste('disp',dispTable(cds)[cond],sep='_')]
    counts <- counts(cds)[, col]
    sizeFactors <- sizeFactors(cds)[col]
    mus <- rowMeans(t(t(counts)/sizeFactors))
    Varofmu <- mus* sum(1/sizeFactors)+disps*mus^2
    Varofmu <- Varofmu/sum(col)^2
    Sigofmu <- sqrt(Varofmu)
    rpkm <- mus / length_bases * 10^3
    Sigofrpkm <- Sigofmu /length_bases *10^3
    out<-data.frame(rpkm=rpkm,Sigofrpkm=Sigofrpkm)
    colnames(out)<-c('rpkm','rpkmstdev')
    return(out)
}


###################################################
### code chunk number 17: s_functions.Rnw:313-339
###################################################
f.barplotwitherror <- function(height, dheight.neg = NULL,
                               dheight.plus = dheight.neg, horiz=FALSE, ...){
    #heith, dheight.neg, dheight.plus are data.frames or coerceable to matrices
    #Each row corresponds to say, a gene, and each column, say, to a condition.
    # The output will be a bar plot with black error bars height-neg, height+pos
    barpositions <- barplot(t(as.matrix(height)), horiz=horiz, beside=TRUE, ...)
    if (!(is.null(dheight.neg) && is.null(dheight.plus)) ){
        if (horiz==FALSE) {
            f.drawerrbar <- function (barpos, h, dh_n, dh_p) {
                segments(x0=barpos, y0= h-dh_n, y1= h+dh_p)
            }
        }
        else {
            f.drawerrbar <- function (barpos, h, dh_n, dh_p) {
                segments(y0=barpos, x0= h-dh_n, x1= h+dh_p)
            }
        }
        if (is.null(dheight.neg)) dheight.neg <- array(0,dim(dheight.plus))
        if (is.null(dheight.plus)) dheight.plus <- array(0,dim(dheight.neg))
        invisible( mapply(f.drawerrbar, barpositions,
                          t(as.matrix(height)),
                          t(as.matrix( dheight.neg)),
                          t(as.matrix( dheight.plus)))
                  )
    }
}


###################################################
### code chunk number 18: s_functions.Rnw:348-358
###################################################
f.getFromPreprocBed <- function( desiredlist,  patternname, samplelist,
                                preprocdir ='./intermediate_results/'){
    covdata <- list()
    for (samplename in samplelist) {
        datafilefullpath <- paste0(preprocdir, patternname, '_', samplename)
        tempfile <- subsetProcessedFiles( desiredlist, datafilefullpath )
        covdata[[samplename]] <- getDataFromCovTemp(tempfile)
    }
    return(covdata)
}


###################################################
### code chunk number 19: s_functions.Rnw:364-380
###################################################
f.plotselectedCovs <- function( covdata, x = seq.int(-10000,9999,10),
                               xlab = 'bp-TSS'){
    par(mfrow=c( length(covdata), length(rownames(covdata[[1]]))))
    par(mar = c(2, 4, 2, 1))
    for (s in seq(along = names(covdata))) {
        sname <- names(covdata)[s]
        for (i in seq(along = rownames(covdata[[sname]]))) {
            id <- rownames(covdata[[sname]])[i]
            plot( x, covdata[[sname]][id,], type='l',
                 main = ifelse( s == 1, id, ''),
                 xlab = ifelse( s == length(covdata), xlab, ''),
                 ylab = ifelse( i == 1, sname, '')
                 )
        }
    }
}


###################################################
### code chunk number 20: s_functions.Rnw:388-403
###################################################
f.avgFromPreprocBed <- function( gs_list,  patternname, samplelist,
                                preprocdir ='./intermediate_results/'){
    x <- list()
    for (samplename in samplelist) {
        message(paste('Processing', samplename))
        datafilefullpath <- paste0(preprocdir, patternname, '_', samplename)
        covdata <- sumCoverageBatch( gs_list, datafilefullpath )
        for (i in seq(along=gs_list)) {
            x[[samplename]][[i]] <-  covdata$sumcovs[[i]] /
                length( covdata$genessummed[[i]] )
        }
        names(x[[samplename]]) <- names(gs_list)
    }
    return(x)
}


###################################################
### code chunk number 21: s_functions.Rnw:406-422
###################################################
f.avgFromPreprocBed_old <- function( gs_list,  patternname, samplelist,
                                preprocdir ='./intermediate_results/'){
    x <- list()
    for (samplename in samplelist) {
        message(paste('Processing', samplename))
        x[[samplename]]<-list()
        for (s in seq(along = gs_list)) {
            datafilefullpath <- paste0(preprocdir, patternname, '_', samplename)
            tempfile <- subsetProcessedFiles( gs_list[[s]], datafilefullpath )
            x[[samplename]][[names(gs_list)[s]]] <-
                sumCoverageFromTemp(tempfile)/length(gs_list[[s]])
            unlink( tempfile )
        }
    }
    return(x)
}


###################################################
### code chunk number 22: s_functions.Rnw:428-450
###################################################
library(MASS)
library(caTools)
f.chipControlWeights <- function( gene_intensity , genesubset, plotresults = FALSE,
                                 windowtoavg = 51) {
    gene_rank <- rank(gene_intensity, ties.method = "random" )
    fitparams <- fitdistr(gene_intensity[genesubset], 'log-normal')
    fitparams <- fitparams$estimate
    if (plotresults) {
        plot.ecdf( log10(gene_intensity[genesubset]), xlim=c(-3,3),
                  xlab = 'log10 Intensity',ylab='CDF',main='chipControl')
        lines( seq(-3,3,0.1), plnorm(10^seq(-3,3,0.1),
                                     fitparams[1],fitparams[2]),lty=3)
    }
    fittedcdf <- plnorm(gene_intensity, fitparams[1], fitparams[2])
    gene_weight <- gene_intensity # Initialize. Will replace values soon.
    gene_weight[ order(gene_rank) ] <- diff( c( 0, fittedcdf[ order(gene_rank)]) )
    # Smooth the weights with a moving window average
    gene_weight[ order(gene_rank) ] <- caTools::runmean(
                                                gene_weight[ order(gene_rank) ],
                                                k= windowtoavg, endrule = 'mean')
    return(gene_weight)
}


###################################################
### code chunk number 23: s_functions.Rnw:455-471
###################################################
f.controlFromPreprocBed <- function( gs_list, gene_intensity,
                                    patternname, samplelist,
                                    preprocdir ='./intermediate_results/'){
    # Function assumes that gene_intensity has names (the gene names)
    genenames <- names(gene_intensity)
    x <- list()
    for (samplename in samplelist) {
        message(paste('Processing', samplename))
        datafilefullpath <- paste0(preprocdir, patternname, '_', samplename)
        weights_list <- lapply(gs_list, function (x)
                               f.chipControlWeights( gene_intensity, x))
        x[[samplename]] <-
            weightedAvgCovBatch( datafilefullpath, weights_list, genenames)
    }
    return(x)
}


###################################################
### code chunk number 24: s_functions.Rnw:474-493
###################################################
f.controlFromPreprocBed_old <- function( gs_list, gene_intensity,
                                    patternname, samplelist,
                                    preprocdir ='./intermediate_results/'){
    # Function assumes that gene_intensity has names (the gene names)
    genenames <- names(gene_intensity)
    x <- list()
    for (samplename in samplelist) {
        x[[samplename]]<-list()
        for (s in seq(along = gs_list)) {
            message(paste('Processing', samplename,names(gs_list)[s]))
            datafilefullpath <- paste0(preprocdir, patternname, '_', samplename)
            tempfile <- subsetProcessedFiles( genenames, datafilefullpath )
            weights <- f.chipControlWeights( gene_intensity, gs_list[[s]] )
            x[[samplename]][[names(gs_list)[s]]] <-
                weightedAvgCovFromTemp( tempfile, weights, genenames)
        }
    }
    return(x)
}


###################################################
### code chunk number 25: s_functions.Rnw:499-519
###################################################
f.plotavgCovs <- function( covdata, covdatacontrol, ylims = NULL,
                          x = seq.int(-10000,9999,10), xlab = 'bp-TSS'){
    par(mfrow=c( length(covdata), length(names(covdata[[1]])) ) )
    par(mar = c(2, 4, 2, 1))
    for (s in seq(along = names(covdata))) {
        sname <- names(covdata)[s]
        for (i in seq(along = names(covdata[[sname]])) ) {
            id <- names(covdata[[sname]])[i]
            plot( x, covdata[[sname]][[id]], type='l',
                 main = ifelse( s == 1, id, ''),
                 xlab = ifelse( s == length(covdata), xlab, ''),
                 ylab = ifelse( i == 1, sname, ''),
                 ylim = c( 0, ifelse( is.null( ylims[[sname]]),
                 1.5*max(covdata[[sname]][[id]]), ylims[[sname]])),
                 col ='red'
                 )
            lines( x, covdatacontrol[[sname]][[id]], type='l', lty=2, col='black')
        }
    }
}


###################################################
### code chunk number 26: s_functions.Rnw:525-534
###################################################
f.weightedcdf <- function( weights, Yvalues ) {
	# arguments must be vectors with names. 
	commonnames <- intersect(names(weights),names(Yvalues))
	w <- weights[commonnames]
	Y <- Yvalues[commonnames] 
	w <- w/sum(w);
	sortorder <- order(Y);
	return( data.frame( val = Y[sortorder], cdf.atval = cumsum( w[sortorder])))
}


###################################################
### code chunk number 27: s_functions.Rnw:543-553
###################################################
f.plotTFcdfs <- function(tf, tf_datacolumns = colnames(tf), genes.in.red, genes.in.blue, ...){
    for (tfname in tf_datacolumns) {
        plot.ecdf(tf[,tfname], col = 'black', do.points = FALSE,
                  main = tfname, ylab = 'cdf', ...)
        plot.ecdf(tf[genes.in.red , tfname], col ='red',
                  do.points = FALSE, add = TRUE)
        plot.ecdf(tf[genes.in.blue , tfname], col = 'blue',
                  do.points = FALSE, add = TRUE)
    }
}


###################################################
### code chunk number 28: s_functions.Rnw:574-588
###################################################
f.fixedcolsplit <- function(string, pattern, names){
    require(stringr)
    outputstr <- list()
    for (i in seq(along = names)){
    	matchpos <- regexpr( pattern = pattern, string, fixed=TRUE )
	outputstr[[i]] <- ifelse( matchpos != -1,
	    substr(string, start = 1, stop = matchpos -1), string)
	string <- ifelse( matchpos != -1, 
	    str_sub(string, start = matchpos + 1), '')
    }
    names(outputstr) <- names
    outputstr[['stringsAsFactors']] = FALSE
    do.call( data.frame, outputstr)
}


###################################################
### code chunk number 29: s_functions.Rnw:594-603
###################################################
f.colsplitandmerge <- function( x, coltosplit, pattern, names) {
    strtosplit <- x[,coltosplit]
    splitcols <- f.fixedcolsplit( strtosplit, pattern, names)
    outdf <- x
    for (nm in names){
	outdf[[nm]] <- splitcols[[nm]]
    }
    return(outdf)
}


###################################################
### code chunk number 30: s_functions.Rnw:606-610
###################################################
f.previewmelted <- function( melteddata){ 
    lapply( split( melteddata[,c('meta','variable')], 
	    melteddata$source, drop = TRUE), table)	
}


###################################################
### code chunk number 31: s_functions.Rnw:613-624
###################################################
f.subsetmelted <- function( m.data, subsetlist) {
    rowstokeep.overall <- rep(FALSE, length(m.data[,1]))
    for (slist in subsetlist){
	rowstokeep <- rep(TRUE, length(rowstokeep.overall))
	for (varname in names(slist)){
	    rowstokeep <- rowstokeep & (m.data[,varname] %in% slist[[varname]])
	}
	rowstokeep.overall <- rowstokeep.overall | rowstokeep
    }
    return( m.data[rowstokeep.overall,] )
}


###################################################
### code chunk number 32: s_functions.Rnw:627-641
###################################################
f.quantXcat_melted <- function( quantdata, catdata) {
    outdat <- dlply( catdata, 'meta', function(x) 
    	subset(quantdata, quantdata[,'GeneSymbol'] %in% x$value))
    idcols <- names(quantdata)
    idcols <- idcols[ idcols != 'value' ]
    outdf <-melt( outdat, id = idcols)
    if (is.factor(catdata$value)) {
	outdf$GeneSymbol <- factor( outdf$GeneSymbol, levels(catdata$value))
    }
    if (is.factor(catdata$meta)) {
	outdf$L1 <- factor(outdf$L1, levels(catdata$meta))
    }
    return(outdf)
}


###################################################
### code chunk number 33: s_functions.Rnw:644-648
###################################################
f.addunitbox <- function() {
    boxdata <- data.frame( xs = c(1,1,-1,-1,1), ys = c(1,-1,-1,1,1))
    geom_path( data=boxdata, aes(x=xs,y=ys,fill=NULL,shape=NULL),color='black')
}


###################################################
### code chunk number 34: s_functions.Rnw:654-659
###################################################
f.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
        paste(toupper(substring(s, 1,1)), substring(s, 2),
	          sep="", collapse=" ")
		  }


