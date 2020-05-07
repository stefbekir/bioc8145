In this tutorial, we will step through using `MuTect2` for somatic SNV and indel calling. `MuTect2` combines the DREAM challenge-winning somatic genotyping engine of the original `MuTect` with the haplotype-based `GATK HaplotypeCaller` backend. 

## 0. Setup

We will use [GATK](https://gatk.broadinstitute.org/hc/en-us) to work with this genomic data. It is available on `Rivanna` as a module, so let us load it.

```bash
module load gatk/4.1.6.0 
```

Let's create a directory where we are going to run all our analyses for today.

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

