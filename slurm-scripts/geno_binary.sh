#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J binary
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --mem=8G

module load R

Rscript scripts/geno_binary.R $SLURM_ARRAY_TASK_ID
