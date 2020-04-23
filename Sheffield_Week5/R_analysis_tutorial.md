# R analysis of ATAC-seq peaks basic tutorial

## Introduction and setup

On rivanna we can load up R like this:

```
module load gcc/8.3.0
module load intel/18.0 intelmpi/18.0 R/3.6.0
```

Next, navigate to a folder where you want to run this tutorial:
```
cd atacseq_tutorial
```

We need 2 pieces of information. First, download this set of differential ATAC-seq peaks published at GEO accession [GSE148396](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148396). You can grab it directly from the command line like this:

```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148396/suppl/GSE148396%5FEBNA2%5Fdependent%5Fopen%5Fchromatin%5Ffiltered%2Ebed%2Egz -O EBNA2_peaks.bed.gz
```

Check out the data to make sure it's what we expect:

```
zcat EBNA2_peaks.bed.gz | head 
```

How many peaks are in this file?

```
zcat EBNA2_peaks.bed.gz | wc -l
```

Notice, we can't use `wc -l` directly on the zipped file because it doesn't have lines. We have to decompress it first, which we are doing with `zcat`.

*Answer question 7*

To load this into R with `rtracklayer` it will be convenient to trim the file down to 3 columns. This type of modification will also be useful for putting this file into GREAT. Do it like this:

```
zcat EBNA2_peaks.bed.gz | cut -f1,2,3 > EBNA2_peaks_3col.bed
```

We're also going to explore the tissue specificity, so you can download this pre-built matrix with tissue specificity:

```
wget http://big.databio.org/open_chromatin_matrix/openSignalMatrix_hg19_quantileNormalized_round4.txt.gz
```

## Start the R analysis

Now let's fire up R:

```
R
```

We'll be using the GenomicDistributions package, which is a package under development in my lab to help us visualize some properties of these regions. There are lots of other R packages that do similar types of analysis, and you could follow their vignettes, but I find this one the easiest to use. Install it with:

```
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
devtools::install_github("databio/GenomicDistributions")
```

If you don't have `devtools` installed, you may need to install it first with `install.packages("devtools")`.



Now, let's load up the package and read in our query regions:

```
library("GenomicDistributions")
queryFile = "EBNA2_peaks_3col.bed"
query = rtracklayer::import(queryFile)
query
```

First a gut check -- how many peaks are there in this R object. 

```
length(query)
```

Now let's see the sizes of these regions

```
summary(width(query))
```

Strange -- there are some really narrow peaks in there! Let's take a closer look at that. We can plot the widths as a histogram.

```
w = width(query)
plotQTHist(w)
```

Here's what you would get out of base R:

```
hist(w)
```

*Answer question 8*


## Chromsome partition plot

Let's first make a chromsome partition plot:
```
x = calcChromBinsRef(query, "hg19")
plotChromBins(x)
```

No red flags there, these peaks are distributed pretty evenly across the chromosomes. Next, we'll take a look at how they distribute relative to annotated features:

```
gp = calcPartitionsRef(query, "hg19")
plotPartitions(gp)
ep = calcExpectedPartitionsRef(query, "hg19")
plotExpectedPartitions(ep)
```




## TSS distance

A common measure of ATAC-seq quality is to make sure the peaks map mostly to TSSs. Here, this is just a subset of the peaks, so it isn't really useful as a quality control, but may still be interesting to look at:

```
TSSdist = calcFeatureDistRefTSS(query, "hg19")
plotFeatureDist(TSSdist, featureName="TSS")
```

## Tissue specificity

Next, we'll look into the tissue specificity of these differential peaks. Here, we have to first load in the cell type specificity matrix file we downloaded earlier.

```
exampleCellMatrixFile = "openSignalMatrix_hg19_quantileNormalized_round4.txt.gz"
cellMatrix = data.table::fread(exampleCellMatrixFile)
```

Then, we can use the GenomicDistributions functions to calculate and plot the tissue specificity:

```
op = calcOpenSignal(query, cellMatrix)
plotOpenSignal(op)
```

*Answer question 9*

If you want to save this to a file, you can do that like this:

```
pdf("my_plot.pdf")
plotOpenSignal(op)
dev.off()
```

*Answer question 10*
