
gencode.all = read.table('gencode.v33.annotation.gtf', sep='\t', header =F)

# the following few ines take some time to execute.
# this could be quicker with the use of data.table as opposed to data.frame
gencode.gene.names = data.frame(sapply(strsplit(sapply(strsplit(
    as.character(gencode.all[,9]),'gene_name '), "[", 2), ";"), "[", 1))

gencode.id.names = data.frame(sapply(strsplit(sapply(strsplit(
    as.character(gencode.all[,9]),'gene_id '), "[", 2), ";"), "[", 1))

gencode.key = cbind(gencode.gene.names, gencode.id.names)

gencode.key = gencode.key[!duplicated(gencode.key[,2]),]

rownames(gencode.key) = gencode.key[,2]

colnames(gencode.key) = c('gene', 'id')

#optional to save intermediate files
save(gencode.key, file = 'gencode.key.Rdata')

#why do I check the dimensions?
dim(all.counts)
dim(gencode.key)

merged.counts.pre = merge(all.counts, gencode.key, by="row.names", all.x=F)

merged.counts = merged.counts.pre[, c(2:(ncol(merged.counts.pre) - 2))]

rownames(merged.counts) = make.names(merged.counts.pre$gene, unique=TRUE)

#Optional exercise:
#Make this code chunk into a function.
