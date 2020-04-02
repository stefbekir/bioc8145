#!/bin/bash

for i in *rmPCRdup.bam
do
    name=$(echo $i | awk -F".rmPCRdup." '{print $1}')
    echo $name
    htseq-count -r pos -f bam --stranded=yes $i gencode.v33.annotation.gtf > $name.gene.counts.txt
done
