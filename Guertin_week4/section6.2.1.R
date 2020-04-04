
#remove the objects path.dir and file.suffix
rm(path.dir, file.suffix)

#these objects no longer exist
print(path.dir)
print(file.suffix)


#define a function
htseq.to.counts.df <- function(path.dir, file.suffix = '.gene.counts.txt') {
    vec.names = c()
    count = 0
    for (txt in Sys.glob(file.path(path.dir, paste0("*", file.suffix)))) {    
        count = count + 1
        experiment.name = strsplit(strsplit(txt, 
                "/")[[1]][length(strsplit(txt, "/")[[1]])], file.suffix)[[1]][1]
        print(experiment.name)
        if (count == 1) {
            all.counts = data.frame(row.names = read.table(txt)[,1])
        }
        vec.names = c(vec.names, experiment.name)      
        all.counts = cbind(all.counts, data.frame(read.table(txt)[,2]))  
    }    
    all.counts = all.counts[1:(nrow(all.counts) - 5),]
    colnames(all.counts) = vec.names
    return(all.counts)
}

all.counts.from.func = htseq.to.counts.df(path.dir = '/scratch/mjg7y/', 
                              file.suffix = '.gene.counts.txt')

#are the results identical?
identical(all.counts.from.func, all.counts)
