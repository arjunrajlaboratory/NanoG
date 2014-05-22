
library(stringr)
library(plyr)
library(reshape2)



datacolumns <- c('A1','A2','A3','A4','A5',
                 'C1','C2','C3','C4','C5')
conds <- c('ALL','GFP','PLUSS','NEG','DIFF',
           'ALL','GFP','PLUSS','NEG','DIFF')
lincfiles <- paste0('./inputdata/CountslincRNA_', datacolumns) 
lincfiles <- data.frame( sample = datacolumns, fname = lincfiles, conds = conds )
linccounts <- ddply( lincfiles, 'sample', function(x) {
    read.delim( as.character(x$fname), header=FALSE, row.names=NULL, skip=1)
})
linccounts[,5]<-NULL
names(linccounts) <- c('sample','Type','chrloc','Ucount',
    'Length','strand','transcript_id')
linccounts <- subset(linccounts, Type == 'transcript')
linccounts <- cbind(linccounts, 
    colsplit( linccounts$chrloc, '(:|-)', names = c('chr','startpos','endpos')))
linccounts <- dcast( linccounts, 
    'transcript_id + chr + startpos + endpos + Length + strand ~ sample',
    value.var = 'Ucount' )
linccounts <- transform(linccounts, GeneSymbol = str_extract(transcript_id, '^[^\\.]*'))
linccounts <- transform(linccounts, transcriptvar = str_extract(transcript_id, '[^\\.]*$'))

save(file='./intermediate_results/linccounts.RData', list='linccounts')
