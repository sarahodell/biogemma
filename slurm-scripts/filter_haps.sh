#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J corfilter
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a_out-.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a_error.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=4
#SBATCH --mem=33G

module load R

chr=$SLURM_ARRAY_TASK_ID
#chr=10

#echo $SLURM_ARRAY_TASK_ID
#Rscript breakup_haplos.R $SLURM_ARRAY_TASK_ID
Rscript scripts/filter_geno.R $chr
#Rscript breakup_haplos.R $chr
