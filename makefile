VPATH = outputdata intermediate_results inputdata



## Run DESeq testing

nbinomindex.RData : cds.RData | computeDeseq.R

	Rscript computeDeseq.R


## Prepare data for DESeq testing

cds.RData : | prepareForDeseqTest.R fullTable.RData

	Rscript prepareForDeseqTest.R


fullTable.RData : geneTable.RData genenames.RData \
 feature_summary_firstannot prepareFullTable.R

	Rscript prepareFullTable.R


geneTable.RData genenames.RData : | \
 findTranscriptGeneSymbols.R \
 feature_quantifications_A1_L001_transcriptsonly_firstannot

	Rscript findTranscriptGeneSymbols.R



## targets to clear up some of the intermediate results. 

cleanAll : cleanIntermediates cleanOutput


cleanIntermediates : 

	find intermediate_results -follow -type f -print0 | xargs -0 rm


cleanOutput : 

	find outputdata -follow -type f -print0 | xargs -0 rm