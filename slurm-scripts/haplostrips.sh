#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/vgt1
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%j-out.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%j-error.txt
#SBATCH -J haplostrips
#SBATCH -t 12:00:00
#SBATCH --mem 3G
#SBATCH --ntasks=1

module load python

scripts/./haplostrips.sh
