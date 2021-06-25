#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J chisq
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --mem 7G

module load R

chr=$SLURM_ARRAY_TASK_ID

Rscript scripts/selection_scan.R $chr
Rscript scripts/haplo_selection_scan.R $chr
