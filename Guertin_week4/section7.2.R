#we can change the sample conditions/comparison groups
sample.conditions.test = factor(c(rep("Replicates1-3", 3), "Replicate4",
                             rep("Replicates1-3", 3), "Replicate4"), 
                             levels=c("Replicates1-3","Replicate4"))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts,
                as.data.frame(sample.conditions.test), ~ sample.conditions.test)

dds = DESeq(deseq.counts.table)

normalized.counts.rna = counts(dds, normalized=TRUE)

DE.results.rep4 = results(dds)

head(DE.results.rep4)

DE.results.rep4.lattice = 
    categorize.deseq.df(DE.results.rep4, 
                        fdr = 0.1, log2fold = 0.0, 
                        treat = 'Replicate 4')

head(DE.results.rep4.lattice)

activated.rep4 = DE.results.rep4.lattice[DE.results.rep4.lattice$response == 
                                 'Replicate 4 Activated',]
dim(activated.rep4)
dim(DE.results.rep4.lattice[DE.results.rep4.lattice$response == 
                                 'Replicate 4 Repressed',])

# Plot these results function
# I am in the minority in that I prefer lattice to ggplot

ma.plot.lattice <- function(ma.df, filename = 'file.name', 
         title.main = "Differential RNA-seq Expression",
         col = c("grey90",  "grey60", "red" , "blue"))
  {
  pdf(paste("MA_plot_", filename, ".pdf", sep=''), 
      useDingbats = FALSE, width=3.83, height=3.83);
  print(xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
               groups=ma.df$response,
               col= col,
                main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"RNA-seq change"), 
               xlab=expression("log"[10]~"Mean of Normalized Counts"),
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), 
                                 par.ylab.text=list(cex=1.1,font=2))));
  dev.off()
  }

ma.plot.lattice(DE.results.rep4.lattice, filename = 'A549_Replicate_4', 
                title.main = "Differential Expression")

                                        #this next chunk is discontinuous in the vignette

gmt.file = paste0(msigdb.url, 'c2.all.v7.0.symbols.gmt')

sample.conditions.4out = factor(c(rep("Replicates1to3", 3), "Replicate4",
                             rep("Replicates1to3", 3), "Replicate4"), 
                             levels=c("Replicates1to3","Replicate4"))

#same as prior
no.batch.4out = deseq.gsea.workflow(merged.counts,
                  batch = FALSE, sample.conditions.4out, 
                  GO.file = gmt.file, name.analysis = 'Leave Rep 4 out', 
                  treatment = 'Leave Rep 4 Out',
                          fdr = 0.1, log2fold = 0.0)

sample.conditions.3out = factor(c(rep("Replicates124", 2), "Replicate3",
                             rep("Replicates124", 3) , "Replicate3","Replicates124"),
                             levels=c("Replicates124","Replicate3"))

sample.conditions.2out = factor(c("Replicates134", "Replicate2", 
                              rep("Replicates134", 3),"Replicate2",
                              rep("Replicates134",2)))

sample.conditions.1out = factor(c("Replicate1", rep("Replicates2to4", 3),
                              "Replicate1", rep("Replicates2to4", 3)), 
                             levels=c("Replicates2to4","Replicate1"))


no.batch.3out = deseq.gsea.workflow(merged.counts,
                  batch = FALSE, sample.conditions.3out, 
                  GO.file = gmt.file, name.analysis = 'Leave Rep 3 out', 
                  treatment = 'Leave Rep 4 Out',
                          fdr = 0.1, log2fold = 0.0)

no.batch.2out = deseq.gsea.workflow(merged.counts,
                  batch = FALSE, sample.conditions.2out, 
                  GO.file = gmt.file, name.analysis = 'Leave Rep 2 out', 
                  treatment = 'Leave Rep 4 Out',
                          fdr = 0.1, log2fold = 0.0)

no.batch.1out = deseq.gsea.workflow(merged.counts,
                  batch = FALSE, sample.conditions.1out, 
                  GO.file = gmt.file, name.analysis = 'Leave Rep 1 out', 
                  treatment = 'Leave Rep 4 Out',
                          fdr = 0.1, log2fold = 0.0)
