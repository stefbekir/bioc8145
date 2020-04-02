#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH -t 72:00:00
#SBATCH --output=bioc8145_guertin_lecture
#SBATCH --partition=standard
#SBATCH -A bioc8145

module load gcc/7.1.0
module load hisat2/2.1.0
module load bedtools/2.26.0
module load samtools/1.10
module load htslib/1.10.2
module load bioconda/py3.6

./200402_rna_seq.sh
