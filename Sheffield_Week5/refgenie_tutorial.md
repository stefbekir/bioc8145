# Refgenie tutorial

To align with bowtie2, we require a bowtie2 index. There are many ways to get reference genome assembly assets, and you can download them directly from the bowtie2 authors, or from the [Illumina iGenomes project](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

Here, I want to show you an easier way: we'll use [refgenie](http://refgenie.databio.org), a tool developed by my lab.

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

So, let's use that in our alignment command:

```
bowtie2 -p 4 -x $(refgenie seek -c refgenie.yaml hg38/bowtie2_index) -1 fastq/tutorial_r1.fastq.gz -2 fastq/tutorial_r2.fastq.gz -S aligned.sam
```

The advanrage of this is that we don't have to keep track of the path to the index because refgenie manages it for us. An entire lab group could use 1 refgenie config, so that a single copy of all the reference resources can be immediately used by all lab members.

*Answer bonus question (optional)*
