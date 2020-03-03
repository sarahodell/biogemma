#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/plink_files
#SBATCH -J plink
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16G

module load plink

#plink --bfile Biogemma_DHLines_600K_Genotypes_binary --freqx --nonfounders --out ../ld_decay/Biogemma_DHLine_allele_freq

plink --bfile Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs --ld-window-kb 100 --ld-window-r2 0 --ld-window 999999 --out ../ld_decay/Biogemma_DHLines_rsquared











