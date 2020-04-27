Introduction to Single Cell RNA-seq
====
  
## Background readings
10 pts per paper
1. [Current best practices in single-cell RNA-seq analysis: a tutorial](https://www.embopress.org/doi/epdf/10.15252/msb.20188746)
2. [Integrative single-cell analysis](https://www.nature.com/articles/s41576-019-0093-7.pdf)
  
## Starting with single cell RNA-seq analysis
  
Request an allocation (previous assignments have the command for that)

Start by making a directory to run the analysis in. Each person has a prefenrece of how they store data, so do what you prefer, but start thinking in a way that you will be able to keep track of each individual analysis. This is just a suggestion.
  
```
mkdir 2020.04.00.scRNA_seq.tutorial
cd 2020.04.00.scRNA_seq.tutorial
mkdir 00_FASTQS
cd 00_FASTQS
```

Next, download FASTQ files from one of the publicly-available data sets on the 10X Genomics support site. 
This example uses the 1,000 cells data set from mouse neurons, consisting of cortex, hippocampus and subventricular zone of an E18 mouse.

## File download

Download the example fastq files:

```
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar
```

The size of this dataset is 6.6G and takes a few minutes to download.
Since this is a tar file and not a tar.gz file, you don't need the -z argument used in previous tutorials to extract it.

```
tar -xvf neuron_1k_v3_fastqs.tar
cd neuron_1k_v3_fastqs
```

*Question 1:* What is the output? With only this information, how many samples were sequenced and what is the information in each type of FASTQ file (R1, R2, and I1)?
  
** Optional, you can run *FASTQC* in these files **

Next, we need a reference transcriptome to align. For future reference, or if you are going to perform this tutorial in your data (and it is not mouse), 10X has several prebuilt transcriptomes in their [support site](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). To avoid issues with space, I made the mouse genome available at: `/project/bioc8145/week6_Alencar/refdata-cellranger-mm10-3.0.0`

## Alignment
Now, we have everything ready to start alignment.
```
module load cellranger/3.1.0
# In case you want to see the options and parameters that are "easily" editible:
cellranger count --help
```
  
To run cellranger count, you need an `--id`. This can be anything you want. The cellranger will use the `--id` to create an output directory. The `--fastq-path` is, like the name suggesting, the path containing your raw FASTQ files. The `--sample` argument is not necessary in this tutorial (hint to question 1), but it was included in the example. Finally, the last nessary argument needed is the path to the transcriptome reference package. However, we will also add an optional argument in this analysis `--expect-cells=1000`. This last argument is most of the time not used, because you are first letting the software "decide" how many cells you saw.

```
cellranger count --id=run_neuron_1kcells --transcriptome=$GENOMEDIR --fastqs=$DATADIR --sample=neuron_1k_v3 --expect-cells=1000
```
*Question 2:* What are the outputs from the alignment step? Explore the folders to understand what was generated. How many cells were sequenced, what are the `Mean Reads per cell` and the `Median Genes per cell` of this particular experimen?  

**Optional: Download [Loupe Browser](https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser) and you can open the `cloupe.cloupe` file (yes.... that is how they decided to name that file...)

Please, if you have any question regarding any step and/or comments regarding this tutorial, please reach out to me: *gf8kz@virginia.edu* or on the Slack channel
