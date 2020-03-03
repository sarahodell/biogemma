#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/permute_haplo
#SBATCH -J glm_threshold
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 500:00:00
#SBATCH --array=1-50%10
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R

# Permuation 1000 times
pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f2 -d,)"

echo $pheno
echo $env

#if [ $env == "ALL" ]
#then 
#    Rscript find_threshold3.R $pheno $env
#fi

Rscript find_threshold3.R $pheno $env
#Rscript GridLMM_plot_by_chr.R $pheno $env










