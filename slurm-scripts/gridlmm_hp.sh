#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/GridLMM_haplotypes
#SBATCH -J gridlmm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-500%50
#SBATCH --ntasks=4
#SBATCH --mem=15G

module load R

base_list=(blank 10 9 9 10 7 10 6 9 9 10)
#base_list=(blank 6 10 6 7 9 9 9 9 8 7)

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f3 -d,)"

base=$(echo ${base_list[$chr]})
hap_list=$(echo $(seq $base 16))

for hapgrp in ${hap_list[@]}; do
    echo $pheno
    echo $env
    echo $chr
    echo $hapgrp
    if [ $env == "ALL" ]
    then
	Rscript GridLMM_run.R $pheno $chr $hapgrp 4
    else
	Rscript GridLMM_pheno_x_env.R $pheno $env $chr $hapgrp 4
    fi
done









