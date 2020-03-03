#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J genofile
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

export PATH=/home/sodell/bin/bin:$PATH

python wgs_geno2.py biogemma/Biogemma_DHLines_600K_Genotypes.vcf.gz Biogemma_082318/Biogemma_0823_geno










