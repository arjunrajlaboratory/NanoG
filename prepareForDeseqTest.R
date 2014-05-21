


Stangle('s_functions.Rnw');source('s_functions.R')

load('./intermediate_results/fullTable.RData')



library(DESeq)
datacolumns <- c('A1','A2','A3','A4','A5',
                 'C1','C2','C3','C4','C5')
conds <- c('ALL','GFP','PLUSS','NEG','DIFF',
           'ALL','GFP','PLUSS','NEG','DIFF')
cds <- newCountDataSet( fullTable[ , datacolumns], conds)
head( counts(cds), n = 3 )
head( fullTable, n = 3 )



cds <- estimateSizeFactors(cds)

# load('./intermediate_results/readtotalsforRPKM.RData')

# pdf('outputdata/sizeFactorsVsReadTotals.pdf', width = 3, height = 4)
# plot( sizeFactors( cds ), readtotalsforRPKM,
#      xlab = 'DESeq Size Factors', ylab = 'Total reads for RPKM normalization',
#      main = 'Sample size factors',cex = 1)
# dev.off()


# pdf('outputdata/SixGenesNormalizedExpression.pdf')
# par(mfrow=c(3,2))
# genestoplot<-c("Gapdh","Actn4","Actn1","Pou5f1","Sox2","Nanog")
# f.plot_counts_alacarte(vectorofgenes = genestoplot,
#                        countlist = counts(cds, normalized = TRUE),
#                        condlist = conditions(cds),sub = 'size normalized')
# dev.off()




cds <- estimateDispersions( cds , method='per-condition',
    sharingMode='maximum', fitType='parametric')


# The options mean:
# \begin{itemize}
#   \item{\texttt{method='per-condition'}} Samples with the same
#     biological condition only are compared when estimating
#     variance. This is as opposed to a 'blind' method in which variance
#     is estimated from all samples together. That method can be good in
#     situations where there will be only a very small number of hits.
#   \item{\texttt{sharingMode='maximum'}} For every gene, there is an
#     $\alpha$ estimate that comes from comparing its values across
#     replicates. This can be higher or lower than the $\alpha(q)$
#     value predicted by the deseq dispersion fit. The
#     \texttt{'maximum'} sharing mode means that the bigger one of the
#     two is kept as the dispersion estimate for the gene. This is the
#     conservative option.
#   \item{\texttt{fitType='parametric'}} Means that the dispersion fit is
#     of the form $\alpha(q)=\alpha_0+\alpha_1/q$, instead of the local
#     regression from the original deseq paper. This is deseq 's new default.
# \end{itemize}



save(file = './intermediate_results/cds.RData', list = 'cds')





# \subsection{Visualizing the dispersion estimates}

# A helper function similar to that in the deseq vignette for plotting
# sample vs fitted dispersions:

# <<>>=
# f.plotDispEsts_bycond <- function( cds, cond, ... ) {
#     columns <- names( conditions( cds )[ conditions( cds ) == cond])
#     plot( rowMeans( counts( cds, normalized = TRUE )[ , columns]),
#         fitInfo( cds, cond)$perGeneDispEsts, pch ='.', log="xy", ... )
#     xg <- 10^seq( -.5, 5, length.out=300 )
#     lines( xg, fitInfo(cds,cond)$dispFun( xg ), col = "red" )
#     title( main = paste0( 'Dispersion Estimates ',
#                paste0( columns, collapse = ',')))
# }

# f.plotDispEstsonly<-function(cds,conditions_tocompare, ...) {
#     colorsforplot=c('black','red','blue','green','gray')
#     xg <- 10^seq( 0, 5, length.out=300 )
#     i<-1;
#     for (cond in conditions_tocompare) {
#         if (i==1) {
#             plot( xg, fitInfo(cds,cond)$dispFun( xg ),
#                         type='l', log="xy",col=colorsforplot[i], ... )
#         } else {
#             lines( xg, fitInfo(cds,cond)$dispFun( xg ),
#                               col=colorsforplot[i] )
#         }
#         i<-i+1
#     }
#     legend('topright', legend = conditions_tocompare, lty = 1,
#     col = colorsforplot[1:length(conditions_tocompare)])
# }
# @


# The plots in Fig. \ref{fig_disperseplots} show that the fitted
# dispersion curves are actually quite different for the different
# conditions. The dispersion of the GFP- and the GFP++ samples are the
# highest, while the two unsorted populations have the smallest
# dispersions.
# <<disperseplots,include=FALSE,echo=TRUE>>=
# par( mfrow = c( 1, 2 ) )
# f.plotDispEsts_bycond(cds, cond = 'ALL',
#                       xlab = 'mean expression (q)', ylab = 'dispersion')
# f.plotDispEstsonly(cds, conditions_tocompare =
#                  c("ALL","GFP","PLUSS","NEG","DIFF"),
#                    xlab='mean expression (q)',ylab='dispersion')
# @

# \begin{figure}
# \begin{center}
# <<figdispplots,fig=TRUE,echo=FALSE>>=
# <<disperseplots>>
# @
# \end{center}
# \caption{\label{fig_disperseplots} Fitted and raw dispersion
#   ($\alpha$) values}
# \end{figure}


# <<echo=FALSE>>=
# invisible(dev.off())
# @
