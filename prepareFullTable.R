
# Having attempted to assign a gene symbol to every transcript, there
# are still a few that for some reason did not get assigned.

library(xtable)

load('./intermediate_results/genenames.RData')
load('./intermediate_results/geneTable.RData')
load('./intermediate_results/countsTable.RData')

geneTable$Gene_Symbol <- genenames
rownames(geneTable)<-geneTable$first_annot

print(head(geneTable[geneTable$annot_type=='refseq'
               & is.na(geneTable$Gene_Symbol),-1 ],n=3))
print(head(geneTable[geneTable$annot_type=='ucsc'
               & is.na(geneTable$Gene_Symbol),-1 ],n=3))
print(head(geneTable[nchar(as.character(geneTable$Gene_Symbol))==0,-1 ],n=3))


# In the REFSEQ category, most of the genes that are listed as NA occur
# because there are more than one genomic region labeled with the same
# REFSEQ code. The output of the RUM pipeline has a \texttt{[[x]]} where
# \texttt{x=1,2,...}, which is not recognized when you search for gene
# symbol. There is good reason to get rid of these genes. Mainly, they
# have very few uniquely mapped reads, since they are after all duplicated.

# I do not understand yet why the majority of UCSC annotations return a
# \texttt{NA}. Many of those correspond to something if you look them up
# on their genome browser. This may have to do with the version of the
# annotation database I am using. Hopefully I'll get this resolved later.

# The VEGA genes that map to an empty Symbol are not found
# at least in the version of the VEGA database from the biomaRt that I
# got. This is a small minority of VEGA annotations.


geneTable <- geneTable[!is.na(geneTable$Gene_Symbol),]
geneTable <- geneTable[!(nchar(as.character(geneTable$Gene_Symbol))==0),]



# \subsection{Multiple isoforms with the same Symbol}
# Many of the remaining transcripts correspond to the same Gene Symbol,
# mostly as different transcript isoforms.


print(length( geneTable$Gene_Symbol ))
print(length( which( duplicated( geneTable$Gene_Symbol ) ) ))


f.showisoforms <- function( genename, genedf ) {
    xt <- genedf[ genedf$Gene_Symbol %in% genename,
                           c( 'first_annot', 'startpos', 'endpos',
                             'Length', 'Gene_Symbol')]
    print(xt)
}
f.showisoforms('Nanog', geneTable)
f.showisoforms('Klf4', geneTable)
f.showisoforms('Lif', geneTable)

# To inform our decision about how to summarize info for each gene,
# let's see what the expression statistics are like for the multiple
# isoforms. Add a column aggregating all counts.

rownames(countsTable)<-countsTable$first_annot
countsTable<-countsTable[rownames(geneTable),-c(1,2)]
datacolumns <- c("A1","A2","A3", "A4","A5",
                 "C1","C2","C3","C4","C5")
countsTable$allcond.counts <- rowSums(countsTable[,datacolumns])
fullTable<-cbind(geneTable,countsTable)



# Create the by-gene-symbol statistics.

library(data.table)
dt<-data.table(fullTable)
stats.bysymbol<-dt[ , list(r.mean=mean(allcond.counts),
                           r.var=var(allcond.counts),
                           rpb.mean=mean(allcond.counts/Length),
                           rpb.var=var(allcond.counts/Length),
                           Length.mean=mean(Length),
                           Length.var=var(Length),
                           n.entries=length(Length)),
                   by=Gene_Symbol]


pdf('outputdata/isoformReadCountsStatistics.pdf')
par( mfrow = c( 1, 2 ) )
plot( log10( stats.bysymbol$r.mean ),
     sqrt( stats.bysymbol$r.var ) / stats.bysymbol$r.mean,
     pch = '.', xlab = 'log10(mean reads)', ylab = 'std dev/ mean reads',
     main = 'byGeneSymbol\n unnormalized reads',
     sub = 'counts summed over all conditions')
plot( log10( stats.bysymbol$rpb.mean ),
     sqrt( stats.bysymbol$rpb.var ) / stats.bysymbol$rpb.mean,
     pch = '.', xlab = 'log10(mean reads per base)',
     ylab = 'std dev rpb/ mean rpb',
     main = 'byGeneSymbol\n Length-normalized reads',
     sub = 'counts summed over all conditions')
dev.off()

# The figure shows that there is substantial variation in the read
# numbers for the isoforms corresponding to each symbol. Part of that
# variation is due to longer genes having more reads in general,


# Keeping overlapping transcripts form the same gene symbol may create problems for the subsequent
# statistical analysis because their read counts will be correlated.

# When I first did the analysis I used only the longest, then I switched
# momentarily to doing only the shortest. The problem with that is that
# there exist ultra small isoforms listed for some genes, especially in
# the VEGA annotation set. (See some of the examples in the tables
# before). It does not seem reasonable to take those as representative
# of the gene.

# We will therefore proceed in two steps:
# \begin{itemize}
#   \item Prefer Refseq to UCSC to VEGA.
#   \item Within an annotation type keep the longest transcript
# \end{itemize}


dt$annot_type <- factor(x=dt$annot_type, levels=c('vega','ucsc','refseq'),
                        ordered=TRUE)


f.selectisoform <- function( annotation_list, length_list, annot_type) {
    annotation_list[ order(annot_type, length_list, decreasing=TRUE)[1] ]
}

selected_isoforms <- dt[,list(locus_annot=
                              f.selectisoform(first_annot,Length,annot_type)),
                        by=Gene_Symbol]
fullTable<-fullTable[selected_isoforms$locus_annot, ]
rownames(fullTable) <- fullTable$Gene_Symbol

f.showisoforms( c('Nanog', 'Klf4', 'Lif'), fullTable )

save(file='./intermediate_results/fullTable.RData',list='fullTable')

