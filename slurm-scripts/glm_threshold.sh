#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM
#SBATCH -J glm_threshold
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
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

#Haplotype threshold
Rscript GridLMM_haplotypes/permute_haplo/find_threshold.R $pheno $env

#Founder threshold
Rscript GridLMM_founderprobs/permute/find_threshold.R $pheno $env

#600K threshold
Rscript GridLMM_600KSNP/permute/find_threshold.R $pheno $env










