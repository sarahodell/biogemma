#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J snake
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=16
#SBATCH --mem 128G

module load conda3

./submit.sh
