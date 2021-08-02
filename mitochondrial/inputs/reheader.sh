#!/bin/bash
#SBATCH --job-name=reheader
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drkthomp@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=1gb
#SBATCH --time=0:10:00
#SBATCH --account=cschnitzler
#SBATCH --qos=cschnitzler
#SBATCH --output=reheader-%j.log
pwd;hostname;date

module load samtools 

for SAMPLE in "$@";
do

# NC_000886.1 is the name of the mitochondrial contig
samtools reheader header bwa_${SAMPLE}_CheMyd_sorted_MITO.bam > temp
mv temp bwa_${SAMPLE}_CheMyd_sorted_MITO.bam 
rm temp 
date
done
