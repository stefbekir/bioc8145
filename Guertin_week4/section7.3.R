sample.conditions.1to3 = factor(sapply(strsplit(as.character(colnames(merged.counts)), 
              '_rep'), '[', 1)[c(1,2,3,5,6,7)], 
              levels=c("A549_no_treatment","A549_30min100nMdex"))

merged.counts.1to3 = merged.counts[,c(1,2,3,5,6,7)]


deseq.gsea.1to3 = deseq.gsea.workflow(merged.counts.1to3,
                  batch = FALSE, sample.conditions.1to3, 
                  GO.file = 'h.all.v7.0.symbols.plusGRresp.gmt', 
                  name.analysis = 'Dexamethasone', 
                  treatment = 'Dexamethasone',
                          fdr = 0.1, log2fold = 0.0)

ma.plot.lattice(deseq.gsea.1to3[[1]], filename = 'A549_RNA_dex', 
        title.main = "Differential Expression")
