#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J geno_qc
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=4
#SBATCH --mem 4G

module load R

Rscript scripts/genotype_qc.R
