#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/probabilities/geno_probs/images
#SBATCH -J plot
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=2
#SBATCH --mem=15G

module load R

Rscript plot_geno_probs.R $SLURM_ARRAY_TASK_ID
