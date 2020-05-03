The SARS-CoV-2 genome, initially reported on January 12 of this year has been studied extensively with the hope of uncovering useful information about COVID-19. Extensive genomic epidemiology of the virus is available at [Nextstrain](https://nextstrain.org/ncov/global).  The coronavirus is an oily membrane packed with genetic instructions to make millions of copies of itself. The instructions are encoded in 30,000 bases — which the infected cell reads and translates into many kinds of virus proteins. An excellent article on the use of mutations to track the virus was published in [The New York Times](https://www.nytimes.com/interactive/2020/04/30/science/coronavirus-mutations.html). We cannot analyze all the data in class, but we will take a look at using `freebayes` to identify mutations in 6 samples that are publicly available in the Short Read Archive.

### SARS-CoV-2 reference and sequencing datasets

The [SARS-CoV-2 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512) can be downloaded from NCBI. For your convenience, I have downloaded it and made it available on Rivanna at: /project/bioc8145/Ratan\_Week7/sars\_cov\_2/sars\_cov\_2.fa. 

Let's check the first few lines of this genome in FASTA format

```bash
head /project/bioc8145/Ratan_Week7/sars_cov_2/sars_cov_2.fa
```

I have also created the (a) BWA index for this reference genome by running the following sequence of commands, (b) FASTA index to allow efficient access to the reference bases

```bash
cd /project/bioc8145/Ratan_Week7/sars_cov_2/

module load gcc/7.1.0 bwa/0.7.17
bwa index sars_cov_2.fa

module load samtools/1.10
samtools faidx sars_cov_2.fa
```

I downloaded a GFF file of annotations for the genome from [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz) which is available at /project/bioc8145/Ratan\_Week7/sars\_cov\_2/sars\_cov\_2.gff. Let's sort it using `bedtools`, and then index it so that we can easily extract annotations for specific genomic intervals.

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz 

mv GCF_009858895.2_ASM985889v3_genomic.gff.gz sars_cov_2.gff.gz

gzip -d sars_cov_2.gff.gz 

module load bedtools/2.26.0 htslib/1.10.2
bedtools sort -i sars_cov_2.gff > sars_cov_2.sorted.gff
bgzip sars_cov_2.sorted.gff
tabix -p gff sars_cov_2.sorted.gff.gz
```

I have downloaded the FASTQ sequences from six experiments from the Short Read Archive (SRA). 4 of the experiments generated paired-end sequences, and 2 of them generated single-end sequences. Several other datasets are available in other places like [GISAID](https://www.gisaid.org/), but require registration and agreement and/or authorization.


### Assignment

For this assignment, we are going to follow the step in this tutorial:

1. use `BWA-mem`, `samblaster`, and `samtools` to align the data to the reference genome, flag duplicates, and sort the resulting BAM file, as we discussed in class. 
2. look at some statistics and plot the coverage of the genome for the 6 samples. 
3. use `freebayes` to detect variants in these data compared to the reference genome.

Create a folder where you are going to run your analyses

```bash
mkdir sars_cov_2
cd sars_cov_2
```

Let's load the appropriate modules, add a directory to our PATH so we can use access tools that are not available as modules without having to enter the absolute path to the directory, and point R_LIBS so that we can all use the same R packages.

```bash
module load gcc/7.1.0 bwa/0.7.17 R/3.6.2 freebayes/0.9.9 vep/90 samtools/1.10 

export PATH=/project/bioc8145/Ratan_Week7/bin:$PATH

export R_LIBS_USER="/project/bioc8145/Ratan_Week7/r_libs"
```

Now let's align the reads to the reference genome, mark the putative PCR duplicates, sort the resulting BAM, and index it. 

```bash
basedir="/project/bioc8145/Ratan_Week7/sars_cov_2"
samples="SRR11140744 SRR11140746 SRR11140748 SRR11140750 SRR11542288 SRR11542289"

for s in ${samples}; do 
    if [ -f ${basedir}/${s}_2.fastq.gz ]; then
        bwa mem -t 1 -K 100000 -R "@RG\tID:${s}\tSM:${s}" \
            ${basedir}/sars_cov_2.fa \
            ${basedir}/${s}_1.fastq.gz \
            ${basedir}/${s}_2.fastq.gz \
        | samblaster --addMateTags \
        | samtools view -b - \
        | samtools sort -o ${s}.bam -O BAM -

    else
        bwa mem -t 1 -K 100000 -R "@RG\tID:${s}\tSM:${s}" \
            ${basedir}/sars_cov_2.fa \
            ${basedir}/${s}.fastq.gz \
        | samblaster --ignoreUnmated \
        | samtools view -b - \
        | samtools sort -o ${s}.bam -O BAM -
    fi
    samtools index ${s}.bam
done 
```

Let's look at some statistics from one of the BAM files.
```bash
samtools flagstat SRR11140744.bam
```

This will print something like this:
```text
1031117 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
25387 + 0 supplementary
285465 + 0 duplicates
1031117 + 0 mapped (100.00% : N/A)
1005730 + 0 paired in sequencing
502865 + 0 read1
502865 + 0 read2
992432 + 0 properly paired (98.68% : N/A)
1005730 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

This tells us that almost all of the reads from this sample were mapped to the reference genome. Most of them are properly paired, i.e. most paired-end reads align with the expected orientation and within the expected distance of each other. Try looking at the other samples, and see if how they fare. 

The putative PCR duplicate rate for these samples is very high. That is because they are generated by a capture process. Whole-genome sequencing typically results in somewhat uniform coverage. Let's see how the coverage across the genome looks for these samples.

First, let us generate the depth of coverage for all the samples when we include all the reads that aligned to the genome

```bash
files=$(ls -1 *.bam)
samtools depth -a -G UNMAP,SECONDARY,QCFAIL -g DUP ${files} > coverage.including_dups.txt
```

Now let's plot the depth and take a look at the depth of coverage distribution using R.

```R
library(tidyverse)

data <- read_tsv("coverage.including_dups.txt", col_names=F) 
colnames(data) <- c("chrom", "pos", "SRR11140744", "SRR11140746","SRR11140748", "SRR11140750", "SRR11542288", "SRR11542289")

pdf("coverage.pdf")
data %>%
    gather(sample, coverage, -chrom, -pos) %>%
    ggplot(aes(pos, coverage, color=sample, group=sample)) + geom_line()
dev.off()
```

You will have to copy the file `coverage.pdf` to your local machine to open and look at it.

From the plot, it seems that the maximum coverage is 8000. That is suspicious. Let us take a closer look at the options for `samtools` depth by typing the following

```bash
samtools depth -h
```

Do you see the issue? Try increasing the maximum depth allowed by `samtools` and re-plot to see how the depth of coverage looks across the genome.

The coverage characteristics appear to confirm our suspicion that these experiments were done using a capture protocol, and more so, it appears that the same capture protocol has been used for all the samples. 

Let's call the variants using `freebayes`, assuming a ploidy of 1, and assuming a homogeneous mixture of cells. 

```bash
freebayes -f ${basedir}/sars_cov_2.fa -p 1 ${files} > variants.vcf
```

This identifies 5 variants. But when we look at the column labeled 'QUAL' in the VCF file, only 2 variants appear to have a high QUAL value. Remember, a QUAL of 10 means that there is a 1 in 10 chance of an error in our assertion of an alternate allele. Interestingly, for the two variants that have high QUAL values, all samples seem to have the alternate allele. 

When we ran `freebayes` above, we ignored the putative PCR duplicates. But should we? This is debatable. Let us try running the same command, but let us include the PCR duplicates in the calculations. 

```bash
freebayes -f ${basedir}/sars_cov_2.fa -p 1 --use-duplicate-reads ${files} > variants.withdups.vcf
```

Alright, so we still see 2 variants with QUAL > 10.  What other assumptions are we making that we should rethink? What about the assumption that these cells are homogeneous? Let's see how many mutations we see if we assume that we don't know the number of samples in the pool, but we would like to find mutations that are seen in at least 10% of the reads at a loci.

```bash
freebayes -f ${basedir}/sars_cov_2.fa -F 0.1 -C 1 --pooled-continuous  --use-duplicate-reads ${files} > variants.frequency.vcf
```

So, now we see 4 variants with QUAL > 10. In the absence of any external validation, it is difficult to say what the right threshold should be, 10%, 20%, or higher. But when thinking about variant calling, you should always be aware of any assumptions made by the caller. If you have calls using a non-sequencing approach, always use it to model a classifier for your calls. All variant callers create several metrics that can be utilized for this purpose. 