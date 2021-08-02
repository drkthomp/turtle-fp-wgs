#!/bin/bash
#SBATCH --job-name=read_counts
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=96:00:00
#SBATCH --account=cschnitzler
#SBATCH --qos=cschnitzler
#SBATCH --output=read_counts-%j.log
pwd;hostname;date
module load R
Rscript read_counts2.R
date

