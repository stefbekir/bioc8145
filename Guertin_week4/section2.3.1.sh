#!/bin/bash

#no dex rep 1 (as defined by GEO description)
fasterq-dump SRR5093037
#no dex rep2
fasterq-dump SRR5093036
#no dex rep3
fasterq-dump SRR5093038
#no dex rep4
fasterq-dump SRR5093039

#30min 100nM dex rep 1
fasterq-dump SRR5093268
#30min 100nM dex rep 2
fasterq-dump SRR5093265
#30min 100nM dex rep 3
fasterq-dump SRR5093267
#30min 100nM dex rep 4
fasterq-dump SRR5093266




mv SRR5093037_1.fastq A549_no_treatment_rep1_PE1.fastq
mv SRR5093036_1.fastq A549_no_treatment_rep2_PE1.fastq
mv SRR5093038_1.fastq A549_no_treatment_rep3_PE1.fastq
mv SRR5093039_1.fastq A549_no_treatment_rep4_PE1.fastq

mv SRR5093268_1.fastq A549_30min100nMdex_rep1_PE1.fastq
mv SRR5093265_1.fastq A549_30min100nMdex_rep2_PE1.fastq
mv SRR5093267_1.fastq A549_30min100nMdex_rep3_PE1.fastq
mv SRR5093266_1.fastq A549_30min100nMdex_rep4_PE1.fastq

mv SRR5093037_2.fastq A549_no_treatment_rep1_PE2.fastq
mv SRR5093036_2.fastq A549_no_treatment_rep2_PE2.fastq
mv SRR5093038_2.fastq A549_no_treatment_rep3_PE2.fastq
mv SRR5093039_2.fastq A549_no_treatment_rep4_PE2.fastq

mv SRR5093268_2.fastq A549_30min100nMdex_rep1_PE2.fastq
mv SRR5093265_2.fastq A549_30min100nMdex_rep2_PE2.fastq
mv SRR5093267_2.fastq A549_30min100nMdex_rep3_PE2.fastq
mv SRR5093266_2.fastq A549_30min100nMdex_rep4_PE2.fastq

gzip *fastq
