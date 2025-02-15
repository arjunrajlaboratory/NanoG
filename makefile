VPATH = intermediate_results inputdata outputdata


allFigures : figure1  figure2  figureCoherence  figurePLUSS \
 figureReconstitution figurePCA  figureFISH  figureFluidigm


## make figures

figureFluidigm : m.deseq.RData \
 inputdata/fluidigmAvg.csv

	Rscript compareToFluidigm.R

figure1 : makefig1.R  m.deseq.RData \
 sorterA.txt  sorterC.txt

	Rscript makefig1.R 


figure2 : makefig2.R  m.deseq.RData  goannot.RData  cats_hits.RData \
 m.kim.RData  m.rex1.RData  m.hayashi.RData  m.morgani.RData

	Rscript makefig2.R


figureCoherence : makefigCoh.R  m.deseq.RData  goannot.RData

	Rscript makefigCoh.R


figurePLUSS : makefigPLUSS.R  m.deseq.RData

	Rscript makefigPLUSS.R


figureReconstitution : makefigReconstitute.R \
 inputdata/FACSreconstitution.csv

	Rscript makefigReconstitute.R


figureHox : makefigHox.R  m.deseq.RData

	Rscript makefigHox.R


figurePCA : makefigSpca.R  cds.RData

	Rscript makefigSpca.R


figureFISH : makefigfish.R \
 new130727.csv \
 new130729.csv \
 new130821.csv \
 new130903.csv \
 new131008.csv 

	Rscript makefigfish.R


figureLincRNA : makefiglincs.R \
 linccounts.RData \
 updatelincs.RData \
 externaldata/guttman11/supp7.csv \
 m.lincs.RData

	Rscript makefiglincs.R


figureMacArthur : makefigmacarthur.R \
 m.deseq.RData \
 m.mac.RData

	Rscript makefigmacarthur.R


figureTopHits :  m.deseq.RData  m.rpkm.RData  goannot.RData

	Rscript makefigtophits.R


## pertaining to lincRNA

linccounts.RData : prepareLincCounts.R \
 inputdata/CountslincRNA_A1 \
 inputdata/CountslincRNA_A2 \
 inputdata/CountslincRNA_A3 \
 inputdata/CountslincRNA_A4 \
 inputdata/CountslincRNA_A5 \
 inputdata/CountslincRNA_C1 \
 inputdata/CountslincRNA_C2 \
 inputdata/CountslincRNA_C3 \
 inputdata/CountslincRNA_C4 \
 inputdata/CountslincRNA_C5

	Rscript prepareLincCounts.R


nbinomindex.wlincs.RData : | computeDeseqWLincs.R \
 linccounts.RData \
 fullTable.RData

	Rscript computeDeseqWLincs.R
 

m.lincs.RData : prepareMLincs.R | nbinomindex.wlincs.RData

	Rscript prepareMLincs.R


# requires bedtools
updatelincs.RData : prepareUpdateLincs.R \
 linccounts.RData \
 externaldata/guttman11/supp1.csv

	Rscript prepareUpdateLincs.R


## run GO-related tasks

cats_hits.RData : computeGOseq.R  goannot.RData m.deseq.RData

	Rscript computeGOseq.R


goannot.RData : prepareGOAnnot.R  fullTable.RData

	Rscript prepareGOAnnot.R


## extract other expression data form literature

m.rex1.RData : prepareRex1.R  fullTable.RData \
 externaldata/marksetalcell/features_GSM850398 \
 externaldata/marksetalcell/features_GSM850399 \
 externaldata/marksetalcell/transcripts_GSM850398 \
 externaldata/marksetalcell/transcripts_GSM850399

	Rscript prepareRex1.R


m.hayashi.RData : prepareHayashi11.R \
 externaldata/HayashiCell11/GSE30056.txt.labelled \
 externaldata/HayashiCell11/indexGSE30056.txt \
 externaldata/HayashiCell11/GPL1261-14790.txt

	Rscript prepareHayashi11.R



m.morgani.RData : prepareMorgani.R \
 externaldata/Morgani13/indexGSE45182.txt \
 externaldata/Morgani13/GSE45182_RefSeqInformation.txt \
 externaldata/Morgani13/GSM1098619_rpkm_GMEM_HV-_rep1.txt \
 externaldata/Morgani13/GSM1098620_rpkm_GMEM_HV-_rep2.txt \
 externaldata/Morgani13/GSM1098621_rpkm_GMEM_HV+_rep1.txt \
 externaldata/Morgani13/GSM1098622_rpkm_GMEM_HV+_rep2.txt \
 externaldata/Morgani13/GSM1098623_rpkm_2i_HV-_rep1.txt \
 externaldata/Morgani13/GSM1098624_rpkm_2i_HV-_rep2.txt \
 externaldata/Morgani13/GSM1098625_rpkm_2i_HV+_rep1.txt \
 externaldata/Morgani13/GSM1098626_rpkm_2i_HV+_rep2.txt \
 externaldata/Morgani13/GSM1098627_rpkm_2i_LIF-_unsorted_rep1.txt \
 externaldata/Morgani13/GSM1098628_rpkm_2i_LIF-_unsorted_rep2.txt \
 externaldata/Morgani13/GSM1098629_rpkm_2i_LIF+_unsorted_rep1.txt \
 externaldata/Morgani13/GSM1098630_rpkm_2i_LIF+_unsorted_rep2.txt

	Rscript prepareMorgani.R

  
m.mac.RData : prepareMacArthur.R \
 externaldata/MacArthur12/GSE40335.txt.labelled \
 externaldata/MacArthur12/sampledescriptions.txt \
 externaldata/MacArthur12/GPL6246.annot.IDandGeneSymbol

	Rscript prepareMacArthur.R


## extract TF-binding data from literature csv

m.kim.RData : prepareKimEtAl.R  externaldata/KimCell08/mmc2.csv

	Rscript prepareKimEtAl.R


## make a melted (stacked) form of the DESeq test data.

m.deseq.RData : prepareMDeseq.R nbinomindex.RData

	Rscript prepareMDeseq.R



## Run DESeq testing

m.rpkm.RData : prepareMRPKM.R \
 fullTable.RData \
 cds.RData \
 inputdata/readsummaryfeaturequant

	Rscript prepareMRPKM.R

# slow step
nbinomindex.RData : | cds.RData computeDeseq.R

	Rscript computeDeseq.R


cds.RData : prepareForDeseqTest.R \
 | fullTable.RData

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