#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/permute_haplo
#SBATCH -J glm_perm
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=2
#SBATCH --mem=16G

module load R/3.5.2

# Permuation 1000 times
chr=$SLURM_ARRAY_TASK_ID
hap_list=$(echo $(seq 2 16))

pheno="tkw_15"
env="ALL"

#while read env; do

for hapgrp in $hap_list; do
    if [ $env == "ALL" ]
    then
	echo "ALL environments"
	Rscript GridLMM_randomized_all.R $pheno $chr $hapgrp 4 1000
    else
	echo "Specific environment"
	Rscript GridLMM_randomized_pheno_x_env.R $pheno $env $chr $hapgrp 4 1000
    fi
done


#done < env_list.txt

#for hapgrp in $hap_list; do
#    Rscript GridLMM_randomized_pheno_x_env.R $pheno $env $chr $hapgrp 4 1000
#done








