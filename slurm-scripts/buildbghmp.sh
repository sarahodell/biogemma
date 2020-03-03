#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J bghmp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00

python MAGICSim_062018/makehmp3.py












