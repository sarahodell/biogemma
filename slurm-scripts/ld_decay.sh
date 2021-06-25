#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J ld_decay
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=16
#SBATCH --mem=125G

module load R

Rscript scripts/ldheatmap.R $SLURM_ARRAY_TASK_ID
