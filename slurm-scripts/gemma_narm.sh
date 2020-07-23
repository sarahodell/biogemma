#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/probabilities/allele_probs
#SBATCH -J narm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-9%1
#SBATCH --ntasks=8
#SBATCH --mem 62G

chr=$SLURM_ARRAY_TASK_ID

Rscript rm_na.R $chr
