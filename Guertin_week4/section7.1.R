github.week4 ='https://raw.githubusercontent.com/stefbekir/bioc8145/master/Guertin_week4/'

#This should be identical to the workflow that we perfomed outside the function:
#load function deseq.gsea.workflow
source(paste0(github.week4, 'deseq.gsea.workflow.R'))

#take a look at these vectors 
sample.conditions = factor(sapply(strsplit(as.character(colnames(merged.counts)), 
              '_rep'), '[', 1), levels=c("A549_no_treatment","A549_30min100nMdex"))

batch = factor(sapply(strsplit(as.character(colnames(merged.counts)), 
              '_rep'), '[', 2))

gmt.file = 'h.all.v7.0.symbols.plusGRresp.gmt'

#recall merged.counts is the counts table

no.batch.all.four = deseq.gsea.workflow(merged.counts,
                  batch = FALSE, sample.conditions, 
                  GO.file = gmt.file, name.analysis = 'No_batch_all_four_reps', 
                  treatment = 'Dexamethasone',
                          fdr = 0.1, log2fold = 0.0)

head(no.batch.all.four[[1]])
head(no.batch.all.four[[2]])

#Apply a batch effect to each replicate

batch.all = factor(sapply(strsplit(as.character(colnames(merged.counts)),
              '_rep'), '[', 2))

batch.all.four = deseq.gsea.workflow(merged.counts, 
                       batch = batch.all, sample.conditions, 
                       GO.file = gmt.file, name.analysis = 'Four_batches_all_four_reps', 
                       treatment = 'Dexamethasone',
                          fdr = 0.1, log2fold = 0.0)


head(batch.all.four[[1]])
head(batch.all.four[[2]])

#is it reasonable to think that replicate 4 is a single batch and reps 1-3 are another?

batch.test = factor(c('batch1', 'batch1', 'batch1', 'batch2', 
               'batch1', 'batch1', 'batch1', 'batch2'))


batch.all.two = deseq.gsea.workflow(merged.counts, batch = batch.test, sample.conditions, 
                          GO.file = gmt.file, name.analysis = 'Two_batches_all_four_reps', 
                          treatment = 'Dexamethasone', fdr = 0.1, log2fold = 0.0)


head(batch.all.two[[1]])
head(batch.all.two[[2]])
