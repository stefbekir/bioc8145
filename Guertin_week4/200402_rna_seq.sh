#!/bin/bash 

idx=hg38
spliceFile=gencode.v33.splice.txt
processors=10

wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_splice_sites.py
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip hg38.fa.gz
hisat2-build hg38.fa ${idx}

gunzip gencode.v33.annotation.gtf.gz
python3 hisat2_extract_splice_sites.py gencode.v33.annotation.gtf > ${spliceFile}

for i in *PE1.fastq.gz
do
    name=$(echo $i | awk -F"_PE." '{print $1}')
    echo $name
    hisat2 -p ${processors} -x ${idx} --rna-strandness FR --known-splicesite-infile ${spliceFile} \
        -1 ${i} -2 ${name}_PE2.fastq.gz | samtools view -bS - | samtools sort -n - -o $name.sorted.bam
    samtools fixmate -m $name.sorted.bam  - | samtools sort - | samtools markdup -rs - $name.rmPCRdup.bam
    htseq-count -r pos -f bam --stranded=yes $name.rmPCRdup.bam gencode.v33.annotation.gtf > $name.gene.counts.txt
    samtools view -b -f 0x40 $name.rmPCRdup.bam | samtools sort - -o ${name}.PE1.sorted
    samtools view -b -f 0x80 $name.rmPCRdup.bam | samtools sort - -o ${name}.PE2.sorted
    genomeCoverageBed -bg -strand - -split -ibam ${name}.PE1.sorted | sortBed -i - > $name.pe1.minus.bedGraph
    genomeCoverageBed -bg -strand + -split -ibam ${name}.PE1.sorted | sortBed -i - > $name.pe1.plus.bedGraph
    genomeCoverageBed -bg -strand + -split -ibam ${name}.PE2.sorted | sortBed -i - > $name.pe2.minus.bedGraph
    genomeCoverageBed -bg -strand - -split -ibam ${name}.PE2.sorted | sortBed -i - > $name.pe2.plus.bedGraph
    cat $name.pe1.plus.bedGraph $name.pe2.plus.bedGraph | sortBed -i - > $name.merged.plus.bedGraph
    touch temp.txt
    echo "track type=bedGraph name=$name.plus color=255,0,0 visibility=2" >> temp.txt
    cat temp.txt $name.merged.plus.bedGraph > $name.plus.bedGraph
    rm temp.txt      
    cat $name.pe1.minus.bedGraph $name.pe2.minus.bedGraph | sortBed -i - > $name.merged.minus.bedGraph
    touch temp.txt
    echo "track type=bedGraph name=$name.minus color=0,0,255 visibility=2" >> temp.txt
    cat temp.txt $name.merged.minus.bedGraph > $name.minus.bedGraph
    rm temp.txt    
    rm *.merged.*us.bedGraph
    rm *.pe*.*us.bedGraph
    gzip $name.minus.bedGraph
    gzip $name.plus.bedGraph
done
