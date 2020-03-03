#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J buildril
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools

python build_ril3.py RILS_10x10.txt B73xOh43_10samples.vcf c10_markers.txt










