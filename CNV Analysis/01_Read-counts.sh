#!/bin/bash
#SBATCH --job-name=read_counts
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=6:00:00
#SBATCH --account=cschnitzler
#SBATCH --qos=cschnitzler-b
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=drkthomp@ucsc.edu     # Where to send mail
#SBATCH --output=read_counts-%j.log
pwd;hostname;date
module load R
Rscript read_counts2.R
date

