#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J pmap
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
module load bcftools

python pscripts/pmap.py biogemma/BiogemmaFounders.csv qtl2/MAGICSim_062018/MAGICSim_pmap.csv










