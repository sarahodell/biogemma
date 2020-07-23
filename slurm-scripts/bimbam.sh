#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/probabilities/allele_probs
#SBATCH -J bimbam
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10%3
#SBATCH --ntasks=8
#SBATCH --mem 63G

module load R
Rscript ../../../scripts/bimbam.R $SLURM_ARRAY_TASK_ID

