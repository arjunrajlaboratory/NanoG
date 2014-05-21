VPATH = intermediate_results inputdata outputdata




allFigures : figure1


## make figures

figure1 : makefig1.R  m.deseq.RData \
 sorterA.txt  sorterC.txt

	Rscript makefig1.R 


figure2 : makefig2.R  m.deseq.RData  goannot.RData  cats_hits.RData

	Rscript makefig2.R


## run GO-related tasks

cats_hits.RData : computeGOseq.R  goannot.RData m.deseq.RData

	Rscript computeGOseq.R


goannot.RData : prepareGOAnnot.R  fullTable.RData

	Rscript prepareGOAnnot.R



## make a melted (stacked) form of the DESeq test data.

m.deseq.RData : prepareMDeseq.R nbinomindex.RData

	Rscript prepareMDeseq.R



## Run DESeq testing: this is a slow step.

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