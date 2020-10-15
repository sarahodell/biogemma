#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J enrichment
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-5000%100
#SBATCH --ntasks=1
#SBATCH --mem 1G

module load R

rep=$SLURM_ARRAY_TASK_ID
Rscript scripts/enrichments_randomization.R $rep
