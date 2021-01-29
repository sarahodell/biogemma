#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/GridLMM_600KSNP/permute
#SBATCH -J glm_perm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-700
#SBATCH --ntasks=4
#SBATCH --mem=8G

module load R

# Permuation 1000 times
# Permuation 1000 times
pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f3 -d,)"
rep="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f4 -d,)"
echo $pheno
echo $env
echo $chr

if [ $env == "ALL" ]
then
    echo "ALL environments"
    Rscript GridLMM_600K_perm.R $pheno $chr $rep 4 100
else
    echo "Specific environment"
    Rscript GridLMM_600K_pheno_x_env.R $pheno $env $chr $rep 4 100
fi
