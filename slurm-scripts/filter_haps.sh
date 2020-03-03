#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/haplotype_probs
#SBATCH -J corfilter
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=3
#SBATCH --mem=24G

module load R/3.5.2

chr=$SLURM_ARRAY_TASK_ID

echo $SLURM_ARRAY_TASK_ID
#Rscript breakup_haplos.R $SLURM_ARRAY_TASK_ID
Rscript breakup_haplos.R $chr












