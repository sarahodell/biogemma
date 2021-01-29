#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/GridLMM_600KSNP
#SBATCH -J glm_snp
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-70
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f3 -d,)"

echo $pheno
echo $chr

if [ $env == "ALL" ]
then
  echo "ALL environments"
  Rscript GridLMM_600K_all.R $pheno $chr 4
else
  echo $env
  Rscript GridLMM_600K_pheno_x_env.R $pheno $env $chr 4
fi
