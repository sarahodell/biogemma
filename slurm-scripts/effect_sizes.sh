#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/
#SBATCH -J effect_plots
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=2
#SBATCH --mem=15G

module load R

#pheno="male_flowering_d6"
#chr=8

Rscript scripts/effect_sizes_comp_all.R
#while read env; do
#    echo $env
#    Rscript effect_sizes.R $chr $pheno $env
#done < env_list.txt
