#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J pgs
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=8
#SBATCH --mem 23G

module load R

#Rscript scripts/pgs.R ALL 8
#Rscript scripts/founder_pgs.R ALL 8
Rscript scripts/haplo_pgs.R ALL 8
