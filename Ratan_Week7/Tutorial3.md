In this tutorial, we will step through using `MuTect2` for somatic SNV and indel calling. `MuTect2` combines the DREAM challenge-winning somatic genotyping engine of the original `MuTect` with the haplotype-based `GATK HaplotypeCaller` backend. Then we will use `Freebayes` and `DNAcopy` to identify loss-of-heterozygosity events in the tumor sample.

## 0. Setup

We will use [GATK](https://gatk.broadinstitute.org/hc/en-us) and [FreeBayes](https://github.com/ekg/freebayes)  to work with these genomic data. They are available on `Rivanna` as modules, so let us load them.

```bash
module load gatk/4.1.6.0 freebayes/0.9.9
```

I have compiled [vcflib](https://github.com/ekg/vcflib) for you,  and all the binaries are available in the folder `/project/bioc8145/Ratan_Week7/bin`. 

An easy way to use them is to add that directory to the PATH variable. Every time we use an executable, all POSIX compliant systems (this includes  Mac OS and Linux systems), look in all directories pointed to by PATH to search for the binary.

```bash
export PATH=/project/bioc8145/Ratan_Week7/bin:$PATH
```

We are also going to use R to detect LoH events. Lets load R, and make sure that you can use the libraries I have installed.

```bash
module load gcc/7.1.0 R/3.6.2
export R_LIBS_USER=/project/bioc8145/Ratan_Week7/r_libs
```

Now, let's create a directory where we are going to run all our analyses for today.

```bash
mkdir somatic_detection
cd somatic_detection
```

## 1. Somatic SNVs using `MuTect2`

To keep things quick enough for the tutorial, we will use the BAM file of alignments from `chr6` of a paired normal/tumor sample. I have made the BAM file available to everyone at `/project/bioc8145/Ratan_Week7/cell_line` on `Rivanna`. 

I have also downloaded a resource file that contains the location and allele frequency of all identified germline variants on `chr6` from [gnomAD](https://gnomad.broadinstitute.org). This is used as one of the filters in `MuTect2`. This file is available on `Rivanna` here `/project/bioc8145/Ratan_Week7/cell_line/af-only-gnomad.chr6.raw.sites.vcf.gz` for your use.

First, we create a VCF file that contains germline and artifactual sites from the normal samples in our dataset. In our case, we have just one sample, so this is not as effective, but this can help identify and eliminate calling artifacts in cohorts. The alignments are done using `/project/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa` as the reference genome, so we use `-L 6` to restrict all analyses to `chr6`. The default value of `--min-sample-count` option is 2, but we will set to it 1, because we have only one normal sample.


```bash
reference=/project/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
gnomad=/project/bioc8145/Ratan_Week7/cell_line/af-only-gnomad.chr6.raw.sites.vcf.gz
tumor_bam=/project/bioc8145/Ratan_Week7/cell_line/tumor.bam
normal_bam=/project/bioc8145/Ratan_Week7/cell_line/normal.bam


# -max-mnp-distance must be set to zero to avoid a bug in GenomicsDBImport
# uses 4 threads by default and takes ~12 minutes 
gatk Mutect2 -R ${reference} -I ${normal_bam} -max-mnp-distance 0 \
  -O normal.vcf.gz -L 6

# Create a GenomicsDB from the normal Mutect2 calls. -V can be used multiple
# times if you have many normals in your cohort
# scales linearly using ~30 seconds
gatk GenomicsDBImport -R ${reference} --genomicsdb-workspace-path pon_db \
  -V normal.vcf.gz -L 6

# Combine the normal calls using CreateSomaticPanelOfNormals
# takes ~13 seconds
gatk CreateSomaticPanelOfNormals -R ${reference} -V gendb://pon_db \
  -O pon.vcf.gz -L 6 --min-sample-count 1
```

Now that we have the panel of normals, let us run `MuTect2` to detect somatic variants. 

```bash
# multiple tumor and normal samples from the same individual ca be specified
# by use of multiple -I.
# uses 4 threads by default and takes ~16 minutes
gatk Mutect2 -R ${reference} -I ${tumor_bam} -I ${normal_bam} \
     -normal norm --germline-resource ${gnomad} \
     --panel-of-normals pon.vcf.gz -O somatic.vcf.gz -L 6
```

You can also call variant in a tumor where you do not have a paired normal by only specifying `-I ${tumor_bam}` in the above command and removing `-normal norm`. You can also call variants in mitochondria using the `--mitochondria` flag.

After you get the variants from `MuTect2`, you should use `FilterMutectCalls` to filter the raw output. Here, we will just run it using the defaults.

```bash
# takes ~8 seconds
gatk FilterMutectCalls -R ${reference} -V somatic.vcf.gz \
  -O filtered.vcf.gz -L 6
```

The above are the basic steps to be taken in somatic SNV calling. In your project, you will want to use more filters such as those that learn any strand-specific errors in your sequencing protocol or use an estimate of cross-sample contamination.  This GATK Best Practices tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136) explains those filters and how to use them in more detail.

## 2. LoH events using DNAcopy

Loss of heterozygosity (LOH) is a common genetic event in cancer development, and is known to be involved in the somatic loss of wild-type alleles in many inherited cancer syndromes. We will take the following steps to identify LoH events in a tumor.

1. Identify heterozygous variants in the germline
2. Calculate the fraction of reads in the tumor that support each of the alleles at each of those loci.
3. Identify large stretches of the genome where the ratios are similar.

So, first let us use `freebayes` to identify variants in the normal and the tumor genome assuming both of them are diploid.

```bash
freebayes -f ${reference} -r 6 --no-indels --no-mnps --no-complex \
    ${normal_bam} ${tumor_bam} \
| vcffilter -f "QUAL > 20" \
| vcf2tsv -g \
| cut -f 1,2,4,5,48,49,53,56 > snp_calls.tsv
```

Now, in R, we are going to limit ourselves to heterozygous SNPs in the germline. Then we are going to use DNAcopy to segment the genome into regions where the ratio is approximately similar.

```R
library(tidyverse)

data <- read_tsv("snp_calls.tsv") %>% 
        mutate(`#CHROM` = as.character(`#CHROM`),
               AO = as.numeric(AO), RO = as.numeric(RO))

norm <- data %>% filter(SAMPLE == "norm") %>% 
        rename("NAO"="AO", "NGT"="GT", "NRO"="RO") %>% select(-SAMPLE)
tumor <- data %>% filter(SAMPLE == "tumor") %>% 
        rename("TAO"="AO", "TGT"="GT", "TRO"="RO") %>% select(-SAMPLE)

df <- left_join(norm, tumor) %>% 
      filter(NGT == "0/1" & !is.na(TGT)) %>%
      mutate(normal_vaf = NAO/(NAO + NRO), tumor_vaf = TAO/(TAO + TRO)) %>%
      select(-NAO, -NGT, -NRO, -TAO, -TGT, -TRO) %>%
      mutate(normal_vaf = abs(0.5 - normal_vaf),
             tumor_vaf = abs(0.5 - tumor_vaf))

library(DNAcopy)

LOH_CNA.object <- CNA(genomdat = df$tumor_vaf, chrom = df$`#CHROM`, 
                      maploc = df$POS, data.type = 'binary')

# Run segmentation analysis
LOH.segmentation <- segment(LOH_CNA.object)
LOH.segmentation.output <- LOH.segmentation$output

# Extract columns needed to plot segmentation results, change column names
LOH_segments <- select(LOH.segmentation.output, chrom, loc.start, loc.end, seg.mean)
colnames(LOH_segments) <- c("CHROM", "START", "END", "MEAN")

LOH_segments$TOP <- .5 + LOH_segments$MEAN
LOH_segments$BOTTOM <- .5 - LOH_segments$MEAN

pdf("LoH.pdf")
ggplot(LOH_segments, aes(x = START, xend=END, y=MEAN, yend=MEAN)) +
    geom_segment()
dev.off()
```

The segments where the mean is close to 0.5 are region of LoH, i.e these are the regions where we see only one allele.
