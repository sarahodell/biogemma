#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/GridLMM_haplotypes/permute_haplo
#SBATCH -J glm_perm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-500%50
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R

# Permuation 1000 times
pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f3 -d,)"

base_list=(blank 8 7 7 8 6 7 6 8 7 7)
#base_list=(blank 7 10 6 7 9 9 8 9 8 2)
base=$(echo ${base_list[$SLURM_ARRAY_TASK_ID]})
hap_list=$(echo $(seq $base 16))


#pheno="grain_yield_15"
#env="ALL"

echo $pheno
echo $env
echo $chr

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

#pheno="female_flowering_d6"
#env="NERAC_2016_WD"
#chr=$SLURM_ARRAY_TASK_ID

#for hapgrp in $hap_list; do
#    Rscript GridLMM_randomized_pheno_x_env3.R $pheno $env $chr $hapgrp 4 1000
#done
