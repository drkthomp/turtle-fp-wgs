#!/bin/bash
#SBATCH --job-name=mitochondrial_subset
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drkthomp@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=1gb
#SBATCH --time=0:10:00
#SBATCH --account=cschnitzler
#SBATCH --qos=cschnitzler
#SBATCH --output=mitochondrial_subset-%j.log
pwd;hostname;date

BASE=/orange/cschnitzler
DATADIR=${BASE}/NS1865-MMartindale_HFHKVDSXY-Lane1-4/NS1865-MMartindale_HFHKVDSXY-Lane1-4/merged_reads/
REF=${BASE}/turtle_seq_files/2021_GCF_015237465.1_rCheMyd1.pri_genomic.fna

module load samtools 

for SAMPLE in "$@";
do
# NC_000886.1 is the name of the mitochondrial contig
samtools view --write-index --reference ${REF} -o bwa_${SAMPLE}_CheMyd_sorted_MITO.bam ${DATADIR}/Platypus_input_files/bwa_${SAMPLE}_CheMyd_sorted.bam NC_000886.1
date
done
