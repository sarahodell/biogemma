#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J reheader
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00


python reheader.py

    






