#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J get10
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

python get10.py










