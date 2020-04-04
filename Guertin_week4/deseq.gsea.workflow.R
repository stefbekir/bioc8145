
deseq.gsea.workflow <- function(counts.matrix, batch = FALSE, sample.conditions, 
                          GO.file = "h.all.v7.0.symbols.gmt",
                          name.analysis = 'test', treatment = 'Dexamethasone',
                          fdr = 0.1, log2fold = 0.0) {
#not the best if statement because it gives me a warning...     
    if (batch == FALSE) { 
        deseq.counts.table = DESeqDataSetFromMatrix(counts.matrix,
                   as.data.frame(sample.conditions), ~ sample.conditions)
        dds = DESeq(deseq.counts.table)
        rld = rlog(dds, blind=TRUE)
        pca.plot = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
        plotPCAlattice(pca.plot, file = paste0("PCA", name.analysis,'.pdf'))
        DE.results = results(dds)
        DE.results.lattice = 
            categorize.deseq.df(DE.results, 
                        fdr = fdr, log2fold = log2fold, treat = treatment)
        }
#This is the code to incorporate the batch effect into the model:    
    else {
        deseq.counts.table = DESeqDataSetFromMatrix(counts.matrix,
                cbind(as.data.frame(batch), as.data.frame(sample.conditions)), 
                ~ batch + sample.conditions)
        dds = DESeq(deseq.counts.table)
        vsd = vst(dds)
#remove batch effect from PCA        
        assay(vsd) = removeBatchEffect(assay(vsd), vsd$batch)
        pca.plot = plotPCA(vsd, intgroup="sample.conditions", returnData=TRUE)
        plotPCAlattice(pca.plot, file =paste0('PCA_', name.analysis, '.pdf'))
        DE.results = results(dds)
        DE.results.lattice = 
            categorize.deseq.df(DE.results, 
                        fdr = fdr, log2fold = log2fold, treat = treatment)
    }
    gene.list= 2^DE.results.lattice$log2FoldChange * 
                     -log(DE.results.lattice$pvalue, base = 10)
    names(gene.list) = row.names(DE.results.lattice)
#rank order
    gene.list = sort(gene.list, decreasing = TRUE)
#remove NA entries
    gene.list = gene.list[!is.na(gene.list)]
#remove duplicated genes
    gene.list = gene.list[!duplicated(names(gene.list))]
    myGO = gmtPathways(GO.file)
    fgRes = as.data.frame(fgsea(pathways = myGO, 
                           stats = gene.list,
                           minSize=5,
                           maxSize=600,
                           nperm=10000))
    print(head(fgRes[order(fgRes$NES, decreasing =TRUE),]))
#plot most significant two 
    pdf(paste("GSEA_first_", name.analysis, ".pdf", sep=''), 
      useDingbats = FALSE, width=6.83, height=3.83);
    print(
        plotEnrichment(myGO[[fgRes[order(fgRes$NES, decreasing =TRUE),][[1]][1]]],
                       gene.list) + labs(title=fgRes[order(fgRes$NES, decreasing =TRUE),][[1]][1])
    )
    dev.off()
    pdf(paste("GSEA_second_", name.analysis, ".pdf", sep=''), 
      useDingbats = FALSE, width=6.83, height=3.83);
    print(
        plotEnrichment(myGO[[fgRes[order(fgRes$NES, decreasing =TRUE),][[1]][2]]],
                       gene.list) + 
        labs(title=fgRes[order(fgRes$NES, decreasing =TRUE),][[1]][2])
    )
    dev.off()
    print('Significantly activated genes')
    print(nrow(DE.results.lattice[DE.results.lattice$response == 
                                  paste0(treatment, " Activated"),]))
    print('Significantly repressed genes')
    print(nrow(DE.results.lattice[DE.results.lattice$response == 
                                  paste0(treatment, " Repressed"),]))
    return(list(DE.results.lattice, fgRes[order(fgRes$NES, decreasing =TRUE),]))
}
