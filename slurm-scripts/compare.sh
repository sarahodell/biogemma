#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/biogemma
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load R/3.3.1

Rscript compare.R 












