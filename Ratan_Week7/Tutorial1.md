We will step through some basic tasks in alignment and variant calling using a handful of Illumina sequencing data sets. This tutorial is largely inspired by the [tutorial on alignment and variant calling](https://github.com/ekg/alignment-and-variant-calling-tutorial)  developed by Erik Garrison. 

## 0. Setup

We will use several tools to work with this genomic data including:

1. [BWA](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [htslib](https://github.com/samtools/htslib)
4. [freebayes](https://github.com/ekg/freebayes)
5. [vcflib](https://github.com/ekg/vcflib/)
6. [vt](https://github.com/atks/vt.git)
7. [sra-tools](https://github.com/ncbi/sra-tools.git)
8. [samblaster](https://github.com/GregoryFaust/samblaster)


Most of these tools are installed on Rivanna as modules. Let us load them so we can use them throughout this tutorial. 

```bash
module load sratoolkit/2.10.5 samtools/1.10 gcc/7.1.0 bwa/0.7.17 \
            freebayes/0.9.9 htslib/1.10.2 
```

You can now check which modules are loaded 

```bash
module list
```

This should print out the names and versions of the modules loaded and available for you to use. As you must have noticed, we did not include `vcflib`, `samblaster`, and `vt` in the above command. These are not available as modules. In that case, we can compile them ourselves. A lot of bioinformatics tools (including the three we want) can be compiled using the pattern I show below. For example, here is how you can compile `vt`

```bash
git clone --recursive https://github.com/atks/vt.git
cd vt && make
```

The above commands create a binary `vt` which we can now use. I have compiled `vt`, `samblaster`, and `vcflib` for you,  and all the binaries are available in the folder `/project/bioc8145/Ratan_Week7/bin`. 

An easy way to use them is to add that directory to the PATH variable. Every time we use an executable, all POSIX compliant systems (this includes  Mac OS and Linux systems), look in all directories pointed to by PATH to search for the binary.

```bash
export PATH=/project/bioc8145/Ratan_Week7/bin:$PATH
```

Let's create a directory where we will run all our analyses for today.

```bash
mkdir snv_detection
cd snv_detection
```


## 1. Aligning E. Coli data with `bwa mem`

Escherichia coli is a gram-negative bacterium that colonizes the lower gut of animals and survives when released to the natural environment, allowing widespread dissemination to new hosts. Pathogenic E. coli strains are responsible for infections of the enteric, urinary, pulmonary, and nervous systems. E. Coli K12 is a common laboratory strain that has lost its ability to live in the human intestine, is ideal for manipulation, and was one of the earliest organisms to be suggested for whole-genome sequencing. Furthermore, the genome is relatively short, only 4,639,221 base pairs, so it is a great place to start learning about alignment and variant calling.


### E. Coli K12 reference genome

I have downloaded the reference genome for you and it is available at `/project/bioc8145/Ratan_Week7/e_coli/E.coli_K12_MG1655.fa` on Rivanna. You can, of course, download it yourself from NCBI (http://www.ncbi.nlm.nih.gov/nuccore/556503834) too. Remember, this is a circular genome, but we represent it linearly in the FASTA format. That format originates from the FASTA software package (developed by David Lipman and Bill Pearson at UVA) but has now become a universal standard in the field of bioinformatics.  Let's copy it to our working space.

```bash
cp /project/bioc8145/Ratan_Week7/e_coli/E.coli_K12_MG1655.fa .
```



### E. Coli sequencing datasets

For testing alignment, we will get two datasets. The first one is a [recently-submitted sequencing run on a K12 strain from the University of Exeter](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1770413). The second one is a [sequencing data set from a different E. Coli strain](http://www.ncbi.nlm.nih.gov/sra/SRX095630[accn]). This one is [famous for its role in foodborne illness](https://en.wikipedia.org/wiki/Escherichia_coli_O104%3AH4#Infection) and is of medical interest.

We can use the [sratoolkit](https://github.com/ncbi/sratoolkit) to directly pull the sequence data (in FASTQ format) from the archive:

```bash
fasterq-dump SRR1770413
fasterq-dump SRR341549
```

`fasterq-dump` is in the SRA toolkit, and allows direct downloading of data from a particular sequencing run ID. SRA stores data in a particular compressed format that isn't compatible with most downstream tools, so it's necessary to put convert it to [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) for further processing.  For each run, we will get two files (you can get 3 files if some reads are paired and others are single-end), as these are paired-end reads.  The first dataset consists of 2x300 bp Illumina reads, and the second dataset consists of 2x100 bp Illumina reads.

`fasterq-dump` does not zip the resulting FASTQ files, and *ideally* you should do that to save space. But for this tutorial, we are going to leave them as they are. 

### Setting up our reference indexes

#### FASTA file index

First, we'll want to allow tools (such as our variant caller) to quickly access certain regions in the reference. This is done using the FASTA index format, which records the lengths of the various sequences in the reference and their offsets from the beginning of the file.

```bash
samtools faidx E.coli_K12_MG1655.fa
```

Now it's possible to quickly obtain any part of the E. Coli K12 reference sequence. For instance, we can get the 200 bp sequence from position 2000000 to 2000200. We'll use a special format to describe the target region `[chr]:[start]-[end]`.

```bash
samtools faidx E.coli_K12_MG1655.fa NC_000913.3:2000000-2000200
```

This returns a FASTA sequence with the subsequence `NC_000913.3:2000000-2000200`.
```text
>NC_000913.3:2000000-2000200
ATACCAGAAATGGTGTACGGGGTTGAGAAATCGTATTTTTTCTTGCGCTCATCAGAAATG
GTGACCTGATTAATCACCACATCAATACGTTTAGAGTCCAGCGACGCCAGCATACCGTCC
CATTTGGTCGGTTTTAGTGACGCCTCAACGCCAAGATGTTTTGCCAGCTGTTGGGCAAAT
TCCACTTCAAAACCGGTTAAT
```


#### BWA index

BWA uses the [FM-index](https://en.wikipedia.org/wiki/FM-index), which a compressed full-text substring index based around the [Burrows-Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform).

To use this index, we first need to build it:

```bash
bwa index E.coli_K12_MG1655.fa
```

You should see `bwa` generate some information about the build process:

```text
[bwa_index] Pack FASTA... 0.03 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.08 seconds elapse.
[bwa_index] Update BWT... 0.02 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.46 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index E.coli_K12_MG1655.fa
[main] Real time: 1.768 sec; CPU: 1.624 sec
```

And, you should notice index files which have been made using the FASTA file name as prefix:

```bash
ls -rt1 E.coli_K12_MG1655.fa*
```

```text
E.coli_K12_MG1655.fa
E.coli_K12_MG1655.fa.fai
E.coli_K12_MG1655.fa.bwt
E.coli_K12_MG1655.fa.pac
E.coli_K12_MG1655.fa.ann
E.coli_K12_MG1655.fa.amb
E.coli_K12_MG1655.fa.sa
```

`E.coli_K12_MG1655.fa` is the reference we copied, `E.coli_K12_MG1655.fa.fai` is the FASTA index, and all the other files are part of the FM index generated by BWA.

### Aligning the datasets against the reference genome

Here's an outline of the steps we'll follow to align the data from the K12 strain against the K12 reference:

1. use BWA to generate SAM records for each read
2. mark PCR duplicates that result from exact duplication of a template during amplification
3. convert the output to BAM
4. sort the output

Let us try these steps one-by-one

```bash
bwa mem -t 1 -K 10000000 -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa SRR1770413_1.fastq SRR1770413_2.fastq \
    > SRR1770413.sam

samblaster -i SRR1770413.sam > SRR1770413.fldup.sam

samtools view -b SRR1770413.fldup.sam > SRR1770413.fldup.bam

samtools sort -O BAM -o SRR1770413.fldup.sorted.bam SRR1770413.fldup.bam
```

The above approach allows us to generate and inspect the intermediate files, but that is not always necessary unless you want to debug the process. It also creates several files that take up space and can lead to confusion. Instead, we will use [unix pipes](https://en.wikiepdia.org/wiki/Pipeline_%28Unix%29) to stream most of these tools together (see this [nice thread about piping `bwa` and `samtools` together on biostar](https://www.biostars.org/p/43677/) for a discussion of the benefits and possible drawbacks of this).

So let's see how we can leverage pipes for this first dataset

```bash
bwa mem -t 1 -K 10000000 -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa SRR1770413_1.fastq SRR1770413_2.fastq \
  | samblaster \
  | samtools view -b - \
  | samtools sort -O BAM -o SRR1770413.piped.fldup.sorted.bam - 
```

Let's do a quick check on this file. Is it the same as the earlier file. Why do you think this is different?

```bash
md5sum SRR1770413.fldup.sorted.bam

md5sum SRR1770413.piped.fldup.sorted.bam 
```

We can run the same workflow for the next sample. In fact, we can run both of the samples using bash variables to our benefit. Let's see how we will do it

```bash
samples="SRR1770413 SRR341549"

for s in ${samples}; do
  bwa mem -t 1 -K 10000000 -R "@RG\tID:${s}\tSM:${s}" \
    E.coli_K12_MG1655.fa ${s}_1.fastq ${s}_2.fastq \
  | samblaster \
  | samtools view -b - \
  | samtools sort -O BAM -o ${s}.fldup.sorted.bam - 
done
```

Note, that I used `"` instead of `'`. Why? 

It's helpful to add a BAM index to the files. This lets us jump around in them quickly using BAM compatible tools that can read the index. 

```bash
samtools index SRR1770413.fldup.sorted.bam

samtools index SRR341549.fldup.sorted.bam
```

## 2. Calling variants

Now that we have the sorted alignments, we can quickly determine variation against the reference by scanning through them using a variant caller. There are many options, including [samtools mpileup](http://samtools.sourceforge.net/samtools.shtml), [freebayes](https://github.com/ekg/freebayes), and [GATK](https://www.broadinstitute.org/gatk/).

For this tutorial, we'll use `freebayes`. It is a versatile caller that is very easy to set up and run, and it produces a well-annotated VCF output that is suitable for immediate downstream filtering.

### Variant calls with `freebayes`

It's quite easy to use `freebayes` provided you have your BAM file completed. We use `--ploidy 1` to indicate that the sample should be genotyped as haploid.

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 \
	SRR1770413.fldup.sorted.bam > SRR1770413.vcf

freebayes -f E.coli_K12_MG1655.fa --ploidy 1 \
	SRR341549.fldup.sorted.bam > SRR341549.vcf
```

The `SRR1770413` sample has an `A->C` mutation at `NC_000913.3:378700`. What is the status of `SRR341549` at that location? `SRR341549` appears to have something going on at `NC_000913.3:378696` which is close. But the comparison of the two is difficult.

### Joint calling

The ideal way to find differences between them is to call them `jointly`. Calling them jointly can help if we have a population of samples to use to help remove calls from paralogous regions. The Bayesian model in `freebayes` combines the data likelihoods from sequencing data with an estimate of the probability of observing a given set of genotypes under assumptions of neutral evolution and a [panmictic](https://en.wikipedia.org/wiki/Panmixia) population. For instance, [it would be very unusual to find a locus at which all the samples are heterozygous](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle). It also helps improve statistics about observational biases (like strand bias, read placement bias, and allele balance in heterozygotes) by bringing more data into the algorithm.

However, in this context, we only have two samples and the best reason to call them jointly is to make sure we have a genotype for each one at every locus where a non-reference allele passes the caller's thresholds in either sample.

We would run a joint call by dropping in both BAMs on the command line to `freebayes`:

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 \
	SRR1770413.fldup.sorted.bam SRR341549.fldup.sorted.bam > e_colis.vcf
```

This only works if we added the (@RG) flags when we aligned the dataset. This is going to take ~9 min, but now we have a single file where we can now easily see that both samples have an MNP at `NC_000913.3:378696`. We will explore that further when we look at the normalization of variants.

### `bgzip` and `tabix`

We can speed up random access to VCF files by compressing them with `bgzip`, in the [htslib](https://github.com/samtools/htslib) package. `bgzip` is a "block-based GZIP", which compresses files in chunks of lines. This chunking lets us quickly seek to a particular part of the file, and support indexes to do so. The default one to use is `tabix`. It generates indexes of the file with the default name `.tbi`.

```bash
bgzip e_colis.vcf

tabix -p vcf e_colis.vcf.gz
```

Now you can pick up a single part of the file. For instance, we could count the variants in a particular region:

```bash
tabix e_colis.vcf.gz NC_000913.3:1000000-1500000 | wc -l
```

If we want to pipe the output into a tool that reads VCF, we'll need to add the `-h` flag, to output the header as well.

```bash
tabix -h e_colis.vcf.gz NC_000913.3:1000000-1500000 \
| vcffilter -f "DP > 10"
```

The `bgzip` format is very similar to that used in BAM, and the indexing scheme is also similar (blocks of compressed data which we build a chromosome position index on top of).

### Quick stats from VCF files

Several tools and scripts are available to manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. An option we will discuss today is [vcflib](https://github.com/ekg/vcflib/), which provides several utilities to manipulate VCFs, but others include [vcftools](http://vcftools.sourceforge.net/), and [bcftools](http://www.htslib.org/doc/bcftools.html). Here is a quick way to check some statistics from a VCF file.

```
vcfstats e_colis.vcf.gz 
```

Here is the output from this command

```text
total variant sites:	47399
of which 47390 (0.99981) are biallelic and 9 (0.000189877) are multiallelic
total variant alleles:	47408
unique variant alleles:	47599

snps:	40085
mnps:	6903
indels:	611
complex:	146

mismatches:	56329

ts/tv ratio:	2.51529
deamination ratio:	1.01223
biallelic snps:	40028 @ 2.67702

ins/del length frequency distribution
length	ins	del	ins/del
1	203	230	0.882609
2	40	28	1.42857
3	11	21	0.52381
4	8	6	1.33333
5	3	6	0.5
6	6	4	1.5
7	2	1	2
8	2	2	1
9	6	5	1.2
10	2	1	2
11	1	1	1
12	4	2	2
13	2	1	2
14	1	1	1
15		1	
16	1		
17	1	1	1
18	2	1	2
19	1		
20			
21			
22		1	
23			
24			
25			
26	1		
27			
28		1	

insertion alleles / deletion alleles:	0.94586
inserted bases / deleted bases:	1.07346

mnp length frequency distribution
length	count
2	919
3	777
4	3227
5	651
6	182
7	687
8	110
9	41
10	179
11	27
12	13
13	48
14	10
15	7
16	17
17	1
18	
19	4
20	
21	
22	1
23	
24	2
total bases in mnps:	16244
```

## 3. The truth shall set you free

The output from `freebayes`, and for that matter all variant callers output a lot of information for each call that can be used to filter them. Ideally, some external validation information should be used to guide the development of pipelines for processing genomic data. In our case, we don't have a 'truth set', we can only filter on intuition, bulk metrics like QUAL, and with an eye for the particular question, we're interested in. What we want is to know the truth for a particular context, so as to understand if our filtering criteria make sense. For example, if we know that the average coverage in the first sample is 55-fold, and the average coverage in the second sample is 131-fold, we can filter requiring some a variant to have a minimum depth of 20 and a maximum depth of 350.

``` bash
bgzip -dc e_colis.vcf.gz | vcffilter -f "DP > 20 & DP < 350" | vcfstats
```

Total variant sites go down from 47,399 to 47,092.

### The NIST Genome in a Bottle truth set for NA12878

A group at the [National Institute of Standards and Technology](https://en.wikipedia.org/wiki/National_Institute_of_Standards_and_Technology) (NIST) has developed a  truth set based on the [HapMap CEU cell line NA12878](https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=GM12878). It's called the [Genome in a Bottle](https://sites.stanford.edu/abms/giab). In addition to characterizing the genome of this cell line using extensive sequencing and manual curation of inconsistencies found between sequencing protocols, the group actually distributes reference material from the cell line for use in validating sequencing pipelines.

To download the truth set, we can head over to the [Genome in a Bottle FTP site](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/) and pick up the latest release. As of writing this, we're at [GiAB v3.3.2](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/). I have downloaded the highly confident calls and the callable region targets by running the following command. These files are available at /project/bioc8145/Ratan_Week7/human.

### The human reference

For variant calling we will use the reference at `/project/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa` on Rivanna.

### Calling variants in [20p12.1](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A12100001-17900000&hgsid=220600397_Vs2XvVv0rRPE9lPwepHAL4Iq3ndi)

To keep things quick enough for the tutorial, let's grab a little chunk of a NA12878 dataset. Let's use [20p12.1](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A12100001-17900000&hgsid=220600397_Vs2XvVv0rRPE9lPwepHAL4Iq3ndi). We'll use a down-sampled alignment made from Illumina HiSeq sequencing of NA12878 that was used as an input to the [NIST Genome in a Bottle](https://github.com/genome-in-a-bottle) truth set for this sample. (Other relevant data can be found in the [GiAB alignment indexes](https://github.com/genome-in-a-bottle/giab_data_indexes).)

We don't need to download the entire BAM file to do this. `samtools` can download the BAM index (`.bai`) provided it hosted alongside the file on the HTTP/FTP server and then use this to jump to a particular target in the remote file.

```bash
samtools view -b ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam 20:12100000-17900000 > NA12878.20p12.1.30x.bam

samtools index NA12878.20p12.1.30x.bam
```

We can call variants as before. Note that we drop the `--ploidy 1` flag. `freebayes` assumes its input is diploid by default. We can use `bgzip` here (using pipes again) to save the extra command for compression.

```bash
reference="/project/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"

freebayes -f ${reference} NA12878.20p12.1.30x.bam \
| bgzip > NA12878.20p12.1.30x.vcf.gz

tabix -p vcf NA12878.20p12.1.30x.vcf.gz
```

This takes ~1.5 min. to complete.

### Comparing our results to the GiAB truth set

In order to compare, we need to exclude things in our output that are outside the callable region, and then intersect with the truth set. That which we don't see in the truth set, and is also in the callable region should be considered a false positive.

First, we'll prepare a reduced representation of this dataset to match 20p12.1:

```bash
cp /project/bioc8145/Ratan_Week7/human/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_* .

# subset the callable regions to chr20 (makes intersection much faster)
cat HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed | grep ^20 > giab_callable.chr20.bed

# subset the high confidence calls to 20p12.1 and rename the sample to match the BAM
tabix -h HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz 20:12100000-17900000 \
    | sed s/HG001/NA12878/ | bgzip >NIST_NA12878_20p12.1.vcf.gz
tabix -p vcf NIST_NA12878_20p12.1.vcf.gz
```

Now, we can compare our results to the calls to get a list of potentially failed sites.

```bash
vcfintersect -r ${reference} -v -i NIST_NA12878_20p12.1.vcf.gz NA12878.20p12.1.30x.vcf.gz \
    | /project/bioc8145/Ratan_Week7/bin/vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.1.30x.giab_failed.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.giab_failed.vcf.gz
```

We can now examine these using `vcfstats`, or manually by inspecting them either serially:

```bash
zcat NA12878.20p12.1.30x.giab_failed.vcf.gz | less -S
```

... or by looking at loci which fail in `samtools tview`.


### Variant normalization

Let's take a look at the statistics of the file from NIST for this region 

```bash
vcfstats NIST_NA12878_20p12.1.vcf.gz
```

As you can see, it does not a large number of MNPs. Our variant file, on the other hand, has a lot of them. In other words, many of the failed variants are unusual in their normalization. For instance:

```text
20   9575773 .  GAGAG  TATAT  1172.52
```

To ensure that comparisons work correctly, we should "normalize" the variants so that they are represented solely as short indels and SNPs.

There are two main problems:

1. `freebayes` represents short haplotypes in the VCF output
2. indels may not be completely left-aligned, there could be additional bases on the call that should be removed so that it can be represented in normalized form

Finally, the variants in the GiAB set have been normalized using a similar process, and doing so will ensure there are not any discrepancies when we compare.

```bash
vcfallelicprimitives -kg NA12878.20p12.1.30x.vcf.gz \
    | vt normalize -r ${reference} - \
    | bgzip >NA12878.20p12.1.30x.norm.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.norm.vcf.gz
```

Here, `vcfallelicprimitives -kg` decomposes any haplotype calls from `freebayes`, keeping the genotype and site-level annotation. (This isn't done by default because in some contexts doing so is inappropriate.) Then `vt normalize` ensures the variants are left-aligned. This isn't important for the comparison, as `vcfintersect` is haplotype-based, so it isn't affected by small differences in the positioning or description of single alleles, but it is good practice.

We can now compare the results again:

```bash
vcfintersect -r ${reference} -v -i NIST_NA12878_20p12.1.vcf.gz NA12878.20p12.1.30x.norm.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip > NA12878.20p12.1.30x.norm.giab_failed.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.norm.giab_failed.vcf.gz
```

Here, we observe why normalization is important when comparing VCF files. Fortunately, the best package available for comparing variant calls to truth sets, [rtgeval](https://github.com/lh3/rtgeval), addresses exactly this concern and also breaks comparisons into three parts matching the three types of information provided by the VCF file--- positional, allele, and genotype. 

### Hard filtering strategies

The failed list provides a means to examine ways to reduce our false positive rate using post-call filtering. We can look at the failed list to get some idea of what might be going on with the failures.

For example, we can test how many of the failed SNPs are removed by applying a simple quality filter and checking the output file's statistics.

```bash
vcffilter -f "QUAL > 10" NA12878.20p12.1.30x.norm.giab_failed.vcf.gz \
    | vcfstats
```

We might also want to measure our sensitivity from different strategies. To do this, just invert the call to `vcfintersect` by removing the `-v` flag (which tells it to invert):

```bash
vcfintersect -r ${reference} -i NIST_NA12878_20p12.1.vcf.gz NA12878.20p12.1.30x.norm.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip > NA12878.20p12.1.30x.norm.giab_passed.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.norm.giab_passed.vcf.gz
```

Now we can test how many variants remain after using the same filters on both:

```bash
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0" NA12878.20p12.1.30x.norm.giab_passed.vcf.gz | vcfstats
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0" NA12878.20p12.1.30x.norm.giab_failed.vcf.gz | vcfstats
```
