# BIOC8145 Assignment 3: Analysis of ChIP-seq Data #

# Instructions #

**0.** Setup a working directory and download the data
You can create a folder/directory in your scratch storage space to work on this analysis. After you log on to Rivanna:

  	$	cd /scratch/{COMPUTINGID}
  	$	mkdir chipseq_project
  	$	cd chipseq_project

You can confirm by showing your current working directory:

  	$	pwd
  	/scratch/{COMPUTINGID}/chipseq_project

Now let’s download the data files. 

  	$	wget https://faculty.virginia.edu/zanglab/bioc8145/data/SRR3728822.fastq.gz
  	$	wget https://faculty.virginia.edu/zanglab/bioc8145/data/SRR3728865.fastq.gz

There is another copy I stored at my /scratch space that you can also copy over by using cp

 	 $	cp -p /scratch/cz3d/bioc8145/data/SRR3728822.fastq.gz ./
 	 $	cp -p /scratch/cz3d/bioc8145/data/SRR3728865.fastq.gz ./
 	 $	mv SRR3728822.fastq.gz LNCaP_AR_DHT.fastq.gz
 	 $	mv SRR3728865.fastq.gz LNCaP_input.fastq.gz

Then we can rename the data files to reflect their sample information:

 	 $	mv SRR3728822.fastq.gz LNCaP_AR_DHT.fastq.gz
 	 $	mv SRR3728865.fastq.gz LNCaP_input.fastq.gz

If you are interested in looking at what the data look like, you can decompress the file and view it:

 	 $	gunzip LNCaP_AR_DHT.fastq.gz
 	 $	head LNCaP_AR_DHT.fastq

You can count the number of lines of the fastq file. 

  	$	wc -l LNCaP_AR_DHT.fastq

It is divisible by 4. Good. Then for simplicity, we can gzip back to the same format:

 	 $	gzip LNCaP_AR_DHT.fastq

We are good to go! 
For the following analysis, it might be easier to get an immediate feedback from each command you are running, so we will perform the analysis interactively on a node using ijob, instead of submitting each job on to Rivanna using SLURM.

  	$	ijob -A bioc8145 -p standard -c 4 -t 4:00:00
  	$	cd /scratch/{COMPUTINGID}/chipseq_project

This will start an interactive job session and put you to a computing node, with 4 CPUs and a time limit of 4 hours. When you finish your analysis, you can exit the interactive job by type in exit


**1. FastQC**:

  	$	module load fastqc
  	$	fastqc LNCaP_AR_DHT.fastq.gz -o ./
  	$	fastqc LNCaP_input.fastq.gz -o ./

You can download the results to the local computer and view it on a web browser. At your local computer,

  	$	scp rivanna.hpc.virginia.edu:/scratch/cz3d/bioc8145/chipseq_project/*fastqc.* ./

Here is some information about how to interpret fastqc result: [https://www.youtube.com/watch?v=GnWSXwQeJ_U](https://www.youtube.com/watch?v=GnWSXwQeJ_U)


**2. bowtie2** Sequence alignment:
We will first load the necessary modules.

	$	module load gcc
	$	module load bowtie2
	$	module load samtools

After loading the modules, we can run bowtie2 using this command, and produce the output in the SAM format:

	$	bowtie2 -p 4 -x 	/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -U 	LNCaP_AR_DHT.fastq.gz -S LNCaP_AR_DHT.sam
	$	bowtie2 -p 4 -x 	/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -U 	LNCaP_input.fastq.gz -S LNCaP_input.sam

At the end of the job running, you will see some summary information on the screen:

	1000000 reads; of these:
	  1000000 (100.00%) were unpaired; of these:
	    16267 (1.63%) aligned 0 times
	    730127 (73.01%) aligned exactly 1 time
	    253606 (25.36%) aligned >1 times
	98.37% overall alignment rate


**3. macs2** peak calling:
We will call peaks using macs2 (using mostly default parameters) from the bowtie2 output SAM format data:

	$	module load macs2
	$	macs2 callpeak -g hs --bdg -q 0.05 -n AR --outdir ./ -t LNCaP_AR_DHT.sam –c LNCaP_input.sam

You can see the process of the job running, and obtained the output files.

Note: If you encounter any problems for obtaining the read mapping and peak calling results. You can find these sample data files for you to work on the following questions 4-6.

Mapped reads in BED format: https://faculty.virginia.edu/zanglab/bioc8145/data/mapped_reads.bed.gz
macs2 peak calling output data: https://faculty.virginia.edu/zanglab/bioc8145/data/test_peaks.narrowPeak 
				https://faculty.virginia.edu/zanglab/bioc8145/data/test_summits.bed 


**4.** The goal is to discover DNA sequence motifs at the transcription factor binding sites (ChIP-seq peaks). Let’s first find the top 5000 strongest AR binding sites (top 5000 peaks) from the macs2 peak calling result. We can use the unix command sort: 

	$	sort -n -r -k 5 AR_peaks.narrowPeak | head -5000 > AR_peaks_top5k.narrowPeak

-n means we want to sort by numerical values; -r means to sort reversely (largest number on the top); -k 5 means sort by the values in the 5th column. head -5000 means taking the top 5000 lines from the file.

Here are 3 different approaches for motif discovery analysis. You can use any option or try them all if you are interested.

Option 1: **MEME**
Now we will get the DNA sequences from the coordinates of these peaks. Let’s first download the whole genome sequence for human genome version hg38.

	$	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	$	gunzip hg38.fa.gz

You can check the file hg38.fa, which contains DNA sequences of the whole human genome. Next let’s use the getfasta function in bedtools to get the DNA sequences of the peaks:

	$	module load gcc
	$	module load bedtools
	$	bedtools getfasta -fi ../fa/hg38.fa -bed AR_peaks_top5k.narrowPeak > 	AR_peaks_top5k.fa

The output file AR_peaks_top5k.fa will be uploaded to the MEME-ChIP web interface ( http://meme-suite.org/tools/meme-chip ). Note that MEME-ChIP will perform de novo motif discovery. After obtaining the motifs discovered, you can use TomTom ( http://meme-suite.org/tools/tomtom ) to compare the top enriched motifs with known motifs from a database, in order to match the de novo motifs with known TF binding motifs.

Option 2: **HOMER**
HOMER (Hypergeometric Optimization of Motif EnRichment) is another useful tool for motif discovery. Let’s first install the HOMER package through conda, a Python package manager. If you don’t have conda installed, let’s install conda first:

	$	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	$	bash Miniconda3-latest-Linux-x86_64.sh

Then follow the instructions to install conda. Then use conda to install HOMER:

	$	conda install -c bioconda homer

Once HOMER is successfully installed, you can run this command for motif discovery:

	$	findMotifsGenome.pl AR_peaks_top5k.narrowPeak hg38 OUTPUT_DIR/ -size 200

You can find more information here: http://homer.ucsd.edu/homer/ngs/peakMotifs.html 


Option 3: **Cistrome AP**
You can directly perform motif search analysis using a function installed in the Cistrome Analysis Pipeline (http://cistrome.org/ap ). You need to upload the BED format peak data (either peaks narrowPeak or summits.bed data files) to the system. In the left panel, look for 

	CISTROME TOOLBOX -> Integrative Analysis -> MOTIF -> SeqPos motif tool
	Then follow the interactive instructions to perform the analysis.


**5.** We will measure the genome-wide distribution of the macs-identified peaks, i.e. compared the locations of the peaks with the locations of genes, especially the promoter regions of the genes (how close the peaks are to a transcription start site). You can use one of several different methods: 

Option 1: **ChIP-seeker**
ChIP-seeker is an R package on Bioconductor. We will install the packages in R if you don’t have it already. In the R environment, run the installation commands:

	> install.packages("BiocManager")
	> BiocManager::install("ChIPseeker")
	> BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
	> BiocManager::install("clusterProfiler")

Then let’s load the package and the needed genome information, i.e., gene annotation information in the hg38 genome version:

	> library(ChIPseeker)
	> library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	> txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	> library(clusterProfiler)

Then we load the macs2 peaks data:

	> peak <- readPeakFile("AR_peaks.narrowPeak")
	> head(peak)

Then we generate a dataset of promoter regions, defined as +/-3000 bp from TSS, and perform the analysis:

	> promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
	> tagMatrix <- getTagMatrix(peak, windows=promoter)
	> dim(tagMatrix)

We can use a pre-calculated tagMatrix to speed up a little:

	> data("tagMatrixList")
	> tagMatrix <- tagMatrixList[[4]]
	> head(tagMatrix)
	> tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

Option 2: **Cistrome AP**
If you’ve tried to use Cistrome Analysis pipeline (http://cistrome.org/ap ) for Question 4, you could also try another tool on Cistrome AP to get the genomic distribution of the peaks. In the left panel, look for:

	CISTROME TOOLBOX -> Integrative Analysis -> ASSOCIATION STUDY -> CEAS: Enrichment on chromosome and annotation

Follow the instructions in the middle panel to specify the uploaded peak data file as well as parameters (it’s always a good idea to begin with all default parameters), you will get the analysis results.

Option 3: **GREAT**
GREAT (Genomic Regions Enrichment of Annotations Tool)( http://great.stanford.edu/public/html/ ) is another web-based tool to perform this analysis among other functions such as gene ontology analyses. Follow the instructions on the website to upload the peak file (BED format, i.e., either narrowPeak or summits.bed from macs2 output) and specify the genome version (hg38), you will get the analysis report. Although the genomic distribution analysis on GREAT has a default setting of 5kb to TSS instead of 3kb, the general trend is similar.


**6.** In order to calculate the FRiP score, we need to figure out how many sequence reads are in macs identified peak regions. One way to do this is to take advantage of the “intersect” function in bedtools, i.e, identify the entries from one BED file that intersect (overlap) with regions in the other BED file. We first load the modules:

	$	module load gcc
	$	module load samtools
	$	module load bedtools

We will convert the bowtie2 generated mapped read SAM file to BED format that are easier to handle:

	$	samtools view -hSb LNCaP_AR_DHT.sam > LNCaP_AR_DHT.bam
	$	bedtools bamtobed -i LNCaP_AR_DHT.bam > LNCaP_AR_DHT.bed

Then we use intersect function to find all the reads that are overlapped with AR narrow peaks: 

	$	bedtools intersect -a LNCaP_AR_DHT.bed -b AR_peaks.narrowPeak -wa > LNCaP_AR_DHT_inpeak.bed

Finally, count the reads and calculate the FRiP score!

	$	wc -l LNCaP_AR_DHT.bed
	$	wc -l AR_DHT_inpeak.bed
