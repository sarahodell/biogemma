#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318
#SBATCH -J alleles
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=6
#SBATCH --array=1-10%3
#SBATCH --mem=48G

module load R/3.5.1

#echo $SLURM_ARRAY_TASK_ID
Rscript wgs_alleleprobs.R $SLURM_ARRAY_TASK_ID












