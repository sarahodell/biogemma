#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J gemma
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-50%10
#SBATCH --ntasks=3
#SBATCH --mem 23G

module load R

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" gemma/pheno_env_list.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" gemma/pheno_env_list.txt | cut -f2 -d,)"

#pheno="male_flowering_d6"
#env="ALL"

Rscript scripts/gemma_plot.R $pheno $env
