
#you need to order the gene list:
#gene.list= 2^DE.results.all.4reps.lattice$log2FoldChange * 
#                 -log(DE.results.all.4reps.lattice$pvalue, base = 10)

#Modified on 220112 need to check
x =  DE.results.all.4reps.lattice$log2FoldChange 
x[x < 0] <- -1
x[x >= 0] <- 1
gene.list= x *(2^abs(DE.results.all.4reps.lattice$log2FoldChange)) *
                 -log(DE.results.all.4reps.lattice$pvalue, base = 10)


names(gene.list) = row.names(DE.results.all.4reps.lattice)

#rank order
gene.list = sort(gene.list, decreasing = TRUE)

#remove NA entries
gene.list = gene.list[!is.na(gene.list)]

#remove duplicated genes
gene.list = gene.list[!duplicated(names(gene.list))]


head(gene.list)


# what file do you want to use? 
#the more categories you have, the higher the multiple testing corrected test statistic

msigdb.url = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/'

#GMT file options
#GO_file= paste0(msigdb.url, "h.all.v7.0.symbols.gmt")
#GO_file= paste0(msigdb.url, 'msigdb.v7.0.symbols.gmt')
#GO_file= paste0(msigdb.url, 'c5.all.v7.0.symbols.gmt')
#GO_file= paste0(msigdb.url, 'c2.all.v7.0.symbols.gmt')



#positive control 
#check the content (second column) of the file to see where I manually copied it from
system('wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/h.all.v7.0.symbols.gmt')
system('wget https://raw.githubusercontent.com/stefbekir/bioc8145/master/Guertin_week4/GO_positive_ctrl.gmt')
system('cat h.all.v7.0.symbols.gmt GO_positive_ctrl.gmt > h.all.v7.0.symbols.plusGRresp.gmt')

GO_file = 'h.all.v7.0.symbols.plusGRresp.gmt'
myGO = gmtPathways(GO_file)


fgRes = as.data.frame(fgsea(pathways = myGO, 
                           stats = gene.list,
                           minSize=5,
                           maxSize=600,
                           nperm=10000))

fgRes[fgRes$padj < 0.1,][order(fgRes$NES[fgRes$padj < 0.1], decreasing =TRUE),]


plotEnrichment(myGO[["GO_DEX_ACTIVATED_GENES_WITH_GR_IN_A549"]],
      gene.list) + labs(title="GR-bound Dex. Activated Genes in A549 (Reddy, et. al., 2009)")
