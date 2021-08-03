#!/bin/bash
#SBATCH --job-name=cnv_calling
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=4:00:00
#SBATCH --account=cschnitzler
#SBATCH --qos=cschnitzler-b
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=drkthomp@ucsc.edu     # Where to send mail
#SBATCH --output=02_CNV-calling-%j.log
pwd;hostname;date
module load R
Rscript 02_CNV-Calling.R
date

