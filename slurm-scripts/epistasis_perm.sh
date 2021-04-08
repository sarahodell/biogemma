#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J epistasis
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-7000%100
#SBATCH --ntasks=1
#SBATCH --mem 3G

if [ ! -f epistasis/permute/report.txt ];then
  touch epistasis/permute/report.txt;
fi

module load R
#pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" epistasis/permute/env_rep_list.txt | cut -f2 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" epistasis/permute/env_rep_list.txt | cut -f1 -d,)"
rep="$(sed "${SLURM_ARRAY_TASK_ID}q;d" epistasis/permute/env_rep_list.txt | cut -f3 -d,)"

Rscript scripts/epistasis_perm.R $chr $env $rep

echo "permutation $env $chr $rep successful" >> epistasis/permute/report.txt
