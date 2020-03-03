#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/
#SBATCH -J model_comp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 10:00:00
#SBATCH --array=1-50
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R/3.6.0

# Permuation 1000 times
pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f2 -d,)"


Rscript GridLMM_haplo_v_600K.R $pheno $env
Rscript GridLMM_founder_v_haplo.R $pheno $env
Rscript GridLMM_founder_v_600K.R $pheno $env












