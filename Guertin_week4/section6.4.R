
#set the no treatment as the reference level (i.e. what happens upon treatment)
sample.conditions = factor(sapply(strsplit(as.character(colnames(merged.counts)), 
              '_rep'), '[', 1), levels=c("A549_no_treatment","A549_30min100nMdex"))

#Note:
# Although not covered herein,
# I typically use these size factors to normalize bedGraph files.
# After I normalize each file (-scale in genomeCoverageBed), 
# I convert to bigWig using bedGraphToBigWig
# and use bigWigMerge to make a signal representative track per condition. 
# estimateSizeFactorsForMatrix(merged.counts)

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts,
                as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

normalized.counts.rna = counts(dds, normalized=TRUE)

rld = rlog(dds, blind=TRUE)

pca.plot = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(pca.plot, file = 'PCA_A549_dex_lattice.pdf')

#it is good to save the R session as you go. 
#I name the image as a time stamp
save.image(file = '200404_R_mjg_1238pm.Rdata')


#this code is discontinuous in the vignette:

DE.results.all.4reps = results(dds)

head(DE.results.all.4reps)

DE.results.all.4reps.lattice = 
    categorize.deseq.df(DE.results.all.4reps, 
                        fdr = 0.1, log2fold = 0.0, treat = 'Dexamethasone')

head(DE.results.all.4reps.lattice)

activated.all = DE.results.all.4reps.lattice[DE.results.all.4reps.lattice$response == 
                                             'Dexamethasone Activated',]

print(activated.all)
