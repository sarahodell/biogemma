#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/biogemma
#SBATCH -J diploid
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem 32G

python make_diploid.py Biogemma_WGS_Founders_filtered_snps_reheader.vcf.gz Biogemma_WGS_Founders_filtered_snps_diploid.vcf














