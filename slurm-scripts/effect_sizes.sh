#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/GridLMM_founderprobs
#SBATCH -J effect
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=4
#SBATCH --mem=3G

module load R

pheno="male_flowering_d6"
chr=8


while read env; do
    echo $env
    Rscript effect_sizes.R $chr $pheno $env
done < env_list.txt
