#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J epistasis
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --mem 4G

module load R
declare -a env_list=( "ALL" "BLOIS_2014_OPT" "BLOIS_2017_OPT" "GRANEROS_2015_OPT" "SZEGED_2017_OPT" "STPAUL_2017_WD" "NERAC_2016_WD" )
#env=$(echo ${env_list[$SLURM_ARRAY_TASK_ID]})

for env in ${env_list[@]}; do
  echo $env
  Rscript scripts/epistasis_scan.R $SLURM_ARRAY_TASK_ID $env
done
