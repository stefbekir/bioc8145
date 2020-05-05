We will step through some discovery and genotyping of structural variants using LUMPY and SVTYPER. 

## 0. Setup

We will use several tools to work with this genomic data including:

1. [LUMPY](https://github.com/arq5x/lumpy-sv)
2. [SVTYPER](https://github.com/hall-lab/svtyper)
3. [BWA](https://github.com/lh3/bwa)
4. [samtools](https://github.com/samtools/samtools)
5. [htslib](https://github.com/samtools/htslib)
6. [samblaster](https://github.com/GregoryFaust/samblaster)

Some of these tools are installed on Rivanna as modules. Lets load them so the we can use them throughout this tutorial. 

```bash
module load samtools/1.10 gcc/7.1.0 bwa/0.7.17 htslib/1.10.2 
```

You can now check which modules are loaded 

```bash
module list
```

This should print out the names and versions of the modules loaded and available for you to use. I have compiled the latest version of LUMPY and it is available on /project/bioc8145/Ratan_Week7/bin. 

An easy way to use them is to add that directory to the PATH variable. Every time we use an executable, all POSIX compliant systems (this includes  Mac OS and Linux systems), look in all directories pointed to by PATH to search for the binary.

```bash
export PATH=/project/bioc8145/Ratan_Week7/bin:$PATH
```

Besides compiling the executable, I have also updated the configuration file for LUMPY at /project/bioc8145/Ratan_Week7/bin/lumpyexpress.config to use the `samtools`, `samblaster`, and `python` that I want.

In order to get access to SVTYPER, we will need to install it using Anaconda. Anaconda is a free and open-source distribution of the Python and R programming languages for scientific computing. `conda ` is a package manager which is shipped with Anaconda.

```bash
module load anaconda/2019.10-py2.7
conda create --name bioc -c bioconda svtyper pysam numpy
conda activate bioc
```

You might have to run, before you are allowed to activate the 

```bash
conda init bash
```

Now, we have SVTYPER installed. Let's check whether it is actually available to us.

```bash
svtyper -h
```

Lets create a directory where we are going to run all our analyses for today.

```bash
mkdir sv_detection
cd sv_detection
```

## Detect SVs from pre-aligned BAMs

To keep things quick enough for the tutorial, we will use a simulated haploid human dataset, limited to 20p12.1 in the genome. Lets copy it into our current working directory

```
cp /project/bioc8145/Ratan_Week7/simulated/sample.20p12.1.30X.bam* .
mv sample.20p12.1.30X.bam sample.bam
mv sample.20p12.1.30X.bam.bai sample.bam.bai
```

LUMPY integrates discordant read-pairs and split-reads to identify SVs. So let us create a file with the discordant read-pairs. We will use the SAM flags to collect all the discordant pairs. Let us think about what we want. We want to throw away alignments where

1. reads are mapped in proper pair, or
2. read is unmapped, or 
3. the mate is unmapped, or
4. the alignment is secondary, or
5. the read is a putative PCR duplicate.

Lets us go the [Explain SAM flags website](https://broadinstitute.github.io/picard/explain-flags.html), and see which SAM flags we should use.

```bash
samtools view -b -F 1294 sample.bam > sample.discordants.bam
samtools index sample.discordants.bam
```

Lets get the split-reads from this data using a script `extractSplitReads_BwaMem` supplied with LUMPY. It is important to note that this script will only work properly if you have the SAM flag has all the information it is looking for. It requires the SA tag (other canonical alignments in case of chimeric alignments) to be set, MC tag (cigar alignment of the mate), and uses the RNEXT, PNEXT.

```bash
samtools view -h sample.bam \
    | extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.bam
samtools index sample.splitters.bam
```

Now, lets run LUMPY using the provided script `lumpyexpress` which will create a file of variants in the VCF format.

```bash
lumpyexpress \
    -B sample.bam \
    -S sample.splitters.bam \
    -D sample.discordants.bam \
    -o sample.vcf
```

Lets take a look at the variants called by LUMPY on IGV. Download the BAM to your local machines to do this. This will allow you to identify artefacts that you must then figure out a way to filter out. One such way is to genotype the variant calls to see if they show that the likelihood of the observed data is supported by a particular genotype. In order to do that we will use SVTYPER

```bash
svtyper -i sample.vcf -B sample.bam -l sample.bam.json -o sv.gt.vcf
```

## Detect SVs from FASTQ files

Let us assume you were starting from the initial FASTQ files 1.fq and 2.fq. If your reference was reference.fa in FASTA format, you could align using BWA-mem such that several of the above steps can be avoided. Here is how you would align

```bash
bwa mem -t 1 -K 10000000 -R "@RG\tID:foo\tSM:bar" reference.fa 1.fq 2.fq \
  | samblaster -d sample.discordants.sam -s sample.splitters.sam --addMateTags \
  | samtools view -b - \
  | samtools sort -O BAM -o sample.bam

```

Lets set up the files, so they are sorted and ready for LUMPY

```bash
samtools view -b sample.discordants.sam \
| samtools sort -O BAM -o sample.discordants.bam
rm sample.discordants.sam 

samtools view -b sample.splitters.sam \
| samtools sort -O BAM -o sample.splitters.bam
rm sample.splitters.sam
```

And now you are ready to use LUMPY and SVTYPER.


## Detect CNVs using read-depth information

We will now learn how to use read-depth to identify CNVs in this same sample. For this we will use `cn.mops`, primarily due to the ease of use. `cn.mops` is designed for use with multiple samples, but can be used for single sample as well.

Lets load R

```bash
module load gcc/7.1.0 R/3.6.2
```

Let us make sure that you can access the installation of cn.mops that I have

```bash
export R_LIBS_USER=/project/bioc8145/Ratan_Week7/r_libs
```

Now lets run `cn.mops`:

```R
library(cn.mops)

BAMFiles <- c("sample.bam")

bamDataRanges <- getReadCountsFromBAM(BAMFiles, WL=1000)
bamDataRanges

# you do not have to do this step with real data 
bamDataRanges <- bamDataRanges[seqnames(bamDataRanges) == 'chr20']
bamDataRanges <- bamDataRanges[12100:17900,]

res <- singlecn.mops(bamDataRanges)
res

res <- calcIntegerCopyNumbers(res)
res
```

This identifies 9 CNV regions, but again, take a look at the calls in IGV.

