#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J enrichment
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-1000%100
#SBATCH --ntasks=2
#SBATCH --mem 10G

module load R

rep=$SLURM_ARRAY_TASK_ID
Rscript scripts/interld_ft_random.R $rep

#Rscript scripts/enrichments_randomization.R $rep
