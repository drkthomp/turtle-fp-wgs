#!/bin/bash
#SBATCH --job-name=ococo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drkthomp@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --account=cschnitzler
#SBATCH --qos=cschnitzler-b
#SBATCH --output=ococo-%j.log
pwd;hostname;date

module load gcc # required for ococo 
module load ococo 

for SAMPLE in "$@";
do
ococo -q 60 -Q 40 -i bwa_${SAMPLE}_CheMyd_sorted_MITO.bam -f ref.fasta -F ${SAMPLE}.fasta -S ${SAMPLE}.stats -V ${SAMPLE}.vcf 
date
done
