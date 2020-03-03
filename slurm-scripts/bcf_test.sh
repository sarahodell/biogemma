#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J tabix
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load tabix

tabix -p vcf hmpv3_founders/Oh43_hmp321_all.vcf.gz










