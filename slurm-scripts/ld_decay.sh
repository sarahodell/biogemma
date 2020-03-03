#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/
#SBATCH -J ld_decay
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=8
#SBATCH --mem=8G

module load R/3.6.1

Rscript ld_decay.R $SLURM_ARRAY_TASK_ID 8









