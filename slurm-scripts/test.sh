#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10

module load R/3.3.1

Rscript test.R ${SLURM_ARRAY_TASK_ID}












