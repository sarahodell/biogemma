#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/
#SBATCH -J snake
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=16
#SBATCH --mem 128G

module load conda3
source activate base

snakemake -p -j 4

