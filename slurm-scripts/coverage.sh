#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/geno_probs
#SBATCH -J coverage
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --mem 16G
#SBATCH --ntasks=2
#SBATCH --array 1-9%5

module load R/3.5.2


Rscript bg_coverage2.R $SLURM_ARRAY_TASK_ID










