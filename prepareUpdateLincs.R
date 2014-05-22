
bedtoolsSystemCallErrorCode <- system('bedtools --version')
if (bedtoolsSystemCallErrorCode == 127){
    stop("this script needs to execute bedtools. Make sure a directory containing bedtools is on your path.")
}

load('./intermediate_results/linccounts.RData')

mergedlincs <- subset( linccounts, transcriptvar == 'merged',
    select = c('GeneSymbol','chr','startpos','endpos','strand'))
mergedlincs <- mergedlincs[ order(mergedlincs$GeneSymbol,decreasing=FALSE),]

supptable <- read.csv('./externaldata/guttman11/supp1.csv')
supptable <- supptable[, c(1:5)]
supptable <- unique(supptable)

# Overlap the two sets with bedtools and go from there:



suppbed <- cbind(supptable[,c('Chromosome', 'Start', 'End', 'Name')], 
    Score=0, supptable[,c('Strand'),drop=FALSE])
write.table(suppbed, file = 'guttman11supp.bed', quote=FALSE, col.names=FALSE,
    row.names=FALSE, sep = "\t")
mergedbed <- cbind(mergedlincs[,c('chr','startpos','endpos','GeneSymbol')],
    Score = 0, mergedlincs[,'strand',drop=FALSE])
write.table(mergedbed, file = 'guttmannew.bed', quote=FALSE, col.names=FALSE,
    row.names=FALSE, sep = "\t")
commandstr <- paste('bedtools intersect -wa -wb -loj',
    '-a guttman11supp.bed', '-b guttmannew.bed', '>guttmancompare.bed')
system( commandstr )
updatelincs <- read.delim('guttmancompare.bed', header=FALSE)
updatelincs <- updatelincs[, c(4,10)]
names(updatelincs) <- c('oldname','newname')
updatelincs$newname[ updatelincs$newname == '.' ] <- NA


save(file = './intermediate_results/updatelincs.RData', list = 'updatelincs')

unlink(c('guttman11supp.bed', 'guttmancompare.bed', 'guttmannew.bed'))