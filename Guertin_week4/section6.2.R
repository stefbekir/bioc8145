# if on rivanna
#module load gcc openmpi R/3.6.1
#you only need to run the next line once to install the packages
source('install_pkgs.R')

library(lattice)
library(DESeq2)
library(dplyr)
library(fgsea)
library(ggplot2)
library(gage)
library(limma)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

path.dir='/scratch/mjg7y/'
file.suffix='.gene.counts.txt'
setwd(path.dir)

#keep track of the experiment names in this list so I can go back and label the columns
vec.names = c()
#we only want generate a data frame label the rows for the first file we loop through
count = 0
#here we will loop through our suffix-defined files in the defined directory
for (txt in Sys.glob(file.path(path.dir, paste0("*", file.suffix)))) {    
      count = count + 1
#this complicated code simply splits strings to extract the experiment name from the file
      experiment.name = strsplit(strsplit(txt, 
            "/")[[1]][length(strsplit(txt, "/")[[1]])], file.suffix)[[1]][1]
      print(experiment.name)
#only do this for the first file, since count is greater than 1 with every other file      
      if (count == 1) {
#generate a data frame with the gene row names and no columns          
          all.counts = data.frame(row.names = read.table(txt)[,1])
          }
#add the experiment name to a growing list of column names      
      vec.names = c(vec.names, experiment.name)      
#for each file (including the first file) add a column with the counts information      
      all.counts = cbind(all.counts, data.frame(read.table(txt)[,2]))  
}

#the last 5 lines are not gene counts, so delete them
all.counts = all.counts[1:(nrow(all.counts) - 5),]
#name the columns
colnames(all.counts) = vec.names

head(all.counts)

dim(all.counts)
