#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM
#SBATCH -J env_test
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --array=1-7

module load R

declare -a env_list=("blank" "ALL" "BLOIS_2014_OPT" "BLOIS_2017_OPT" "GRANEROS_2015_OPT" "SZEGED_2017_OPT" "STPAUL_2017_WD" "NERAC_2016_WD")
declare -a pheno_list=("male_flowering_d6" "male_flowering_days" "female_flowering_d6" "female_flowering_days" "grain_yield_15" "tkw_15" "total_plant_height" "harvest_grain_moisture")
env=$(echo ${env_list[$SLURM_ARRAY_TASK_ID]})
#pheno=$(echo ${env_list[$SLURM_ARRAY_TASK_ID]})
#echo $env
#echo $pheno

for p in "${pheno_list[@]}"; do
    Rscript environment_test.R $env $p;
done
