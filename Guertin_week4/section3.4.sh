#!/bin/bash

idx=hg38
spliceFile=gencode.v33.splice.txt
processors=3

for i in *PE1.fastq.gz
do
    name=$(echo $i | awk -F"_PE1." '{print $1}')
    echo $name
    hisat2 -p ${processors} -x ${idx} --rna-strandness FR --known-splicesite-infile ${spliceFile} \
        -1 ${i} -2 ${name}_PE2.fastq.gz | samtools view -bS - | samtools sort -n - -o $name.sorted.bam
done
