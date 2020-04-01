#!/bin/bash

for i in *.rmPCRdup.bam
do
    name=$(echo $i |  awk -F".rmPCRdup.bam" '{print $1}')
    echo $name
    samtools view -b -f 0x40 $i | samtools sort - -o ${name}.PE1.sorted
    samtools view -b -f 0x80 $i | samtools sort - -o ${name}.PE2.sorted
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
