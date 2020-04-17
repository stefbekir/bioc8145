# First steps for ATAC-seq analysis

## Introduction and setup

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

*Answer questions 3-1 and 3-2*

Let's run fastqc on our input files
```
cd ..
mkdir fastqc
fastqc --noextract --outdir fastqc fastq/tutorial_r1.fastq.gz
fastqc --noextract --outdir fastqc fastq/tutorial_r2.fastq.gz
```

Now, let's clone the PEPATAC repository, which has some scripts we'll need

```
git clone https://github.com/databio/pepatac.git
```

## Trim adapters (optional)

Next we'll do adapter trimming. This is an optional step, so I won't show you how to install skewer, but you can go do it if you want. The tutorial files have already had their adapters trimmed, but this is how we'd do it if not:

```
skewer -f sanger -t 8 -m pe -x pepatac/tools/NexteraPE-PE.fa --quiet -o fastq/tutorial fastq/tutorial_r1.fastq.gz fastq/tutorial_r2.fastq.gz
```

You may want to run `fastqc` on your input files after trimming, to see how they look and if your trimming was successful.


## Alignment

Make sure you have `bowtie2` in your PATH. You can load it on rivanna like this:

```
module load bowtie2
```

To align with bowtie2, we require a bowtie2 index. There are many ways to get reference genome assembly assets, and you can download them directly from the bowtie2 authors, or from the [Illumina iGenomes project](https://support.illumina.com/sequencing/sequencing_software/igenome.html). Here, I want to show you an easier way: we'll use [refgenie](http://refgenie.databio.org), a tool developed by my lab.

Install refgenie with python how you would typically install a python package:

```
pip install --user refgenie
```

You can read more details in the [refgenie documentation](http://refgenie.databio.org), but for now, we just need to initialize a configuration file and then download the bowtie2 index with `refgenie pull`:

```
refgenie init -c refgenie.yaml
refgenie pull -c refgenie.yaml hg38/bowtie2_index
```

Now we have the index managed by refgenie. We can retrieve the local path to it with:

```
refgenie seek -c refgenie.yaml hg38/bowtie2_index
```

*Answer bonus question (optional)*

So, let's use that in our alignment command:

```
bowtie2 -p 4 -x $(refgenie seek -c refgenie.yaml hg38/bowtie2_index) -1 fastq/tutorial_r1.fastq.gz -2 fastq/tutorial_r2.fastq.gz -S aligned.sam
```

## Check out some alignment statistics

We might be interested in what percentage of our reads aligned to the mitochondria.

```
samtools sort aligned.sam | samtools view -S -b - > aligned.bam
```

```
samtools index aligned.bam
samtools idxstats aligned.bam
samtools idxstats aligned.bam | grep -we 'chrM'
samtools idxstats aligned.bam | grep -we 'chrM' | cut -f 3
samtools idxstats aligned.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd'
```

*Answer question 3-3*


## Shift reads

As we discussed in the lecture, ATAC-seq reads are often shifted to account for a 9bp duplication in the transpose insertion. However, the adjustment is small and for simple peak calling, is probably not that important, and is commonly left out of basic ATAC analysis. Therefore, we'll skip the step here, but if you want to learn how to do it, you could try a tool like [alignmentSieve](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html).

## Peak calling

You can load macs2 on rivanna like this:

```
module load macs2/2.1.2
```

Use the `callpeak` function to call peaks:

```
macs2 callpeak -t aligned.bam -n "tutorial" -g hs -f BAM -q 0.01 --shift 0 --nomodel
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
cut -f1 tutorial_peaks.narrowPeak | uniq -c | sort -k 1 -r
```

*Answer question 3-4*

<!-- For fixed-width peaks, you could use these params: '--shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01' -->
 