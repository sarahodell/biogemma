#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/MAGICSim_121718
#SBATCH -J genofile
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
module load bcftools

python ../impute/qtl2/genofile2.py MAGIC_DHSimAll_121718_chr10.vcf.gz MAGICSim_121718_geno.csv










