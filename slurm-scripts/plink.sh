#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J plink
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=24
#SBATCH --mem 185G

module load plink

#plink --bfile Biogemma_DHLines_600K_Genotypes_binary --freqx --nonfounders --out ../ld_decay/Biogemma_DHLine_allele_freq

plink --threads 24 --bfile genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs inter-chr --out stats/ld_decay/Biogemma_DHLines_rsquared_all_chroms
