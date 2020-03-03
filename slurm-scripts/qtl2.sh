#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/geno_probs
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=4
#SBATCH --mem 32G

module load R/3.5.2

#Rscript qtl2.R 
#Rscript sim_geno.R
#Rscript bg_percentage.R
#Rscript control_file.R
#Rscript get_ibd.R
#Rscript wgs_genoprobs.R
#Rscript recode_ibd.R
#./split_geno.sh
#Rscript bg_coverage.R
#python wgs_founders2.py
#Rscript wgs_founders3.R 10
#Rscript founder_alleles.R $SLURM_ARRAY_TASK_ID

Rscript bg_coverage.R 10

