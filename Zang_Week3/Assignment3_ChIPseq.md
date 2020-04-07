# BIOC8145 Assignment 3

Analysis of ChIP-seq data #

(100 pts total; Due Monday 4/13/2020)

The androgen receptor (AR) is a transcription factor (TF) activated by androgenic hormones (e.g., dihydrotestosterone, DHT) and regulates expression of keys genes in prostate differentiation and function. Aberrant elevation of AR signaling is a crucial driver of prostate cancer. Determining the genome-wide binding pattern of AR is important for understanding the oncogenic gene regulation in prostate cancer cells. In this assignment, you will explore the genomic binding profile of AR in a DHT-treated prostate cancer cell line, LNCaP (Lymph Node Carcinoma of the Prostate), by analyzing some ChIP-seq data using the tools you just learned.

ChIP-seq data to be used in this assignment can be found in a paper published in _Nucleic Acids Research_ (Malinen M et al. 2017. PMID: 27672034). The data can be downloaded from [https://faculty.virginia.edu/zanglab/bioc8145/data/](https://faculty.virginia.edu/zanglab/bioc8145/data/)

You can find 2 data files in this directory:

        SRR3728822.fastq.gz                ChIPseq\_LNCaP\_AR\_DHT (GSM2219854)

        SRR3728865.fastq.gz                ChIPseq\_LNCaP\_input (GSM2219876)

1. First check the sequence read quality of both datasets using **fastqc**. Describe the overall fastqc report and assess the sequence quality. (20 pts)

#Sample commands:

$ module load fastqc

$ fastqc -o _$OUTDIR__SRR3728822.fastq.gz_

2. Now let's map the sequence reads to the human genome (hg38) using **bowtie2**. For each dataset, list the total number of reads, total number of aligned reads, and the overall alignment rate. (30 pts)

Note: bowtie2 index can be found at [https://faculty.virginia.edu/zanglab/bioc8145/bowtie2\_index/](https://faculty.virginia.edu/zanglab/bioc8145/bowtie2_index/)

Copy bowtie2 index to your designated directory (INDEX\_DIR)

#Sample commands:

$ module load gcc

$ module load bowtie2

$ bowtie2 -p 4 -x _${INDEX\_DIR}_/GRCh38 -U _???.fastq.gz_ -S _NAME_.sam

3. Identify AR binding sites (call peaks) from the mapped ChIP-seq reads data using **macs2**. Briefly describe the results. What is the redundant rate of the sequence reads? (Hint: look into the macs2 output data files) How many peaks have you identified? Under what parameter settings? How many peaks have a fold enrichment greater than 5 and 10, respectively? (30 pts)

#Sample commands:

$ module load macs2

$ macs2 callpeak -h

$ macs2 callpeak -t _???.sam_ -c _???.sam_ -n _NAME_

4. Look for DNA sequence motifs enriched in the 5000 strongest AR binding sites using **MEME**. What motifs did you find? Besides AR, did you find other motifs? If you did, what does it indicate? (20 pts)

MEME suite for ChIP-seq analysis can be accessed at[http://meme-suite.org/tools/meme-chip](http://meme-suite.org/tools/meme-chip)

5. (Bonus 1) Can you quantify the genome-wide distribution of the identified peaks? i.e., how many peaks are located in gene promoter regions, say, \&lt; 3kb from any transcription start site (TSS)? How many peaks are located in intronic regions or intergenic regions? Based on the observation, do you think AR is a promoter-binding factor or a distal enhancer-binding factor? (+20 pts)

Hint: You can use an R-package called ChIPseeker [https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html](https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)

6. (Bonus 2) Can you figure out a way to calculate the FRiP score (Fraction of Reads in Peaks) of this AR ChIP-seq dataset? Describe the tool you used. Attached the source code if you wrote your own code to do this. (+40 pts)



Note: Always use SLURM to submit jobs on Rivanna. You should NEVER run a job on the head node.
