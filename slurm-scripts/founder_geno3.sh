#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J genofile
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
module load bcftools

python pscripts/genofile2.py biogemma/BiogemmaFounders_600K_Genotypes.vcf.gz qtl2/MAGICSim_062018/MAGICSim_founder_geno2.csv










