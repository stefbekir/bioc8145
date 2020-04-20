# First steps for ATAC-seq analysis

## Introduction and setup

Request an allocation if you like:

```
ijob -A bioc8145 -p instructional -c 4 -t 1:00:00
```


First we'll create a folder for the tutorial and a subfolder for our raw data:

```
mkdir atacseq_tutorial
cd atacseq_tutorial
mkdir fastq
cd fastq
```

## File download

Download the example fastq files:

```
wget http://big.databio.org/pepatac/tutorial_r1.fastq.gz
wget http://big.databio.org/pepatac/tutorial_r2.fastq.gz
```

Let's take a quick look at what's in these files:

```
zcat tutorial_r1.fastq.gz | head
time zcat tutorial_r1.fastq.gz | head
```

*Answer questions 3 and 4*

Let's run fastqc on our input files
```
cd ..
mkdir fastqc
module load fastqc/0.11.5
fastqc --noextract --outdir fastqc fastq/tutorial_r1.fastq.gz
fastqc --noextract --outdir fastqc fastq/tutorial_r2.fastq.gz
```

Now, let's get the Nextera adapter sequence

```
wget https://raw.githubusercontent.com/databio/pepatac/master/tools/NexteraPE-PE.fa
```

## Trim adapters (optional)

Next we'll do adapter trimming. This is an optional step, so I won't show you how to install skewer, but you can go do it if you want. The tutorial files have already had their adapters trimmed, but this is how we'd do it if not:

I pre-installed [skewer](https://github.com/relipmoc/skewer) for you on rivanna.

```
/project/bioc8145/week5/skewer/skewer -f sanger -t 1 -m pe -x NexteraPE-PE.fa --quiet -o fastq/tutorial fastq/tutorial_r1.fastq.gz fastq/tutorial_r2.fastq.gz
```

You may want to run `fastqc` on your input files after trimming, to see how they look and if your trimming was successful.


## Alignment

Make sure you have `bowtie2` in your PATH. You can load it on rivanna like this:

```
module load gcc/7.1.0
module load bowtie2/2.2.9
```

To align with bowtie2, we require a bowtie2 index. There are many ways to get reference genome assembly assets, and you can download them directly from the bowtie2 authors, or from the [Illumina iGenomes project](https://support.illumina.com/sequencing/sequencing_software/igenome.html). 

If you're on rivanna, I've downloaded the index and placed it at `/project/bioc8145/week5/bowtie2_index`. If you are interested, we recently developed a tool to help in downloading and manage reference indexes called *refgenie*. I posted a [refgenie tutorial](refgenie_tutorial.md) that you may find interesting. 

For now, let's use the existing bowtie2 index in our alignment command:

```
time bowtie2 -p 4 -x /project/bioc8145/week5/bowtie2_index/hg38 -1 fastq/tutorial_r1.fastq.gz -2 fastq/tutorial_r2.fastq.gz -S aligned.sam
```

## Shift reads

As we discussed in the lecture, ATAC-seq reads are often shifted to account for a 9bp duplication in the transpose insertion. However, the adjustment is small and for simple peak calling, is probably not that important, and is commonly left out of basic ATAC analysis. Therefore, we'll skip the step here, but if you want to learn how to do it, you could try a tool like [alignmentSieve](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html).

## Check out some alignment statistics

We might be interested in what percentage of our reads aligned to the mitochondria.

```
module load samtools/1.10
samtools sort aligned.sam | samtools view -S -b - > aligned.bam
```

```
samtools index aligned.bam
samtools idxstats aligned.bam
samtools idxstats aligned.bam | grep -we 'chrM'
samtools idxstats aligned.bam | grep -we 'chrM' | cut -f 3
```

*Answer question 5*

## Peak calling

You can load macs2 on rivanna like this:

```
module load macs2/2.1.2
```

Use the `callpeak` function to call peaks:

```
macs2 callpeak -t aligned.bam -n "tutorial" -g hs -q 0.01 --shift 0 --nomodel
```

Take a quick look at our peaks like this:

```
head tutorial_peaks.narrowPeak
```

How many peaks are there?

```
wc -l tutorial_peaks.narrowPeak
```

How are they distributed across chromosomes?

```
cut -f1 tutorial_peaks.narrowPeak | uniq -c | sort -k 1 -n -r
```

*Answer question 6*

<!-- For fixed-width peaks, you could use these params: '--shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01' -->
 