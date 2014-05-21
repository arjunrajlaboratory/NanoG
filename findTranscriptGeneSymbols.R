
# \subsection{Input data transcript choice}

# The particular transcripts that are in the input data are those chosen
# by the RUM aligner as the transcriptome. This is done by running the
# \verb+rum_indexes+ function and choosing the mm9 mouse build. These seem to
# be an aggregation of all the refseq, ucsc ``known genes,'' and vega
# genes available from the UCSC genome browser.

library(org.Mm.eg.db)
library(biomaRt)


fileforgeneinfo<-"feature_quantifications_A1_L001_transcriptsonly_firstannot"

geneTable<- read.delim(paste0("./inputdata/",fileforgeneinfo),
                       header=TRUE,stringsAsFactors=FALSE)
geneTable<- geneTable[,c(1:7,12)]

save( file = './intermediate_results/geneTable.RData', list = 'geneTable' )


# \subsection{Find gene symbols for each gene}
# Find the gene symbols for each gene. The Refseq and UCSC are found
# from the bioconductor mouse database

# The VEGA annotations require a different procedure. Firstly, the
# database is not included in \texttt{org.Mm.eg.db}, and second, the
# mapping function I have to use for VEGA does not return a
# corresponding row if it doesn't find that particular entry in the
# database. (As opposed to the \texttt{org.Mm.eg.db} methods above,
# which return \texttt{NA} when something isn't found)


compute.getGeneSymbols <- function(geneTable) {
    attach(geneTable)
    genenames<-array("",dim=length(rownames(geneTable)))

    message('Finding REFSEQ Gene Symbols')
    genenames[annot_type=='refseq'] <- select(
              org.Mm.eg.db, keys=first_annot[annot_type=='refseq'],
              cols=c("SYMBOL"),keytype="REFSEQ")$SYMBOL

    message('Finding UCSC Gene Symbols')
    genenames[annot_type=='ucsc']<- select(
              org.Mm.eg.db, keys=first_annot[annot_type=='ucsc'],
              cols=c("SYMBOL"),keytype="UCSCKG")$SYMBOL

    message('Finding VEGA Gene Symbols')
    mart <- useMart('vega')
    dset<-useDataset('mmusculus_gene_vega',mart)
    vegaquery <- getBM(
                 attributes=c('vega_transcript_id','mgi_symbol'),
                 filters='vega_transcript_id',
                 values=first_annot[annot_type=='vega'], mart=dset)
    names(genenames)<-first_annot
    genenames[vegaquery$vega_transcript_id]<-vegaquery$mgi_symbol

    detach("geneTable")

    return(genenames)
}
genenames <- compute.getGeneSymbols(geneTable)

save(file='./intermediate_results/genenames.RData', list = 'genenames')