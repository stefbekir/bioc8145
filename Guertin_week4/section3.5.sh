#!/bin/bash

for i in *.sorted.bam
do
      name=$(echo $i | awk -F".sorted.bam" '{print $1}')
      echo $name
      samtools fixmate -m $i - | samtools sort - | samtools markdup -rs - $name.rmPCRdup.bam
done
