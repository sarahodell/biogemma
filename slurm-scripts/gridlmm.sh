#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM
#SBATCH -J gridlmm
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=15G

module load R

base_list=(blank 6 10 6 7 9 9 9 9 8 7)
base=$(echo ${base_list[$SLURM_ARRAY_TASK_ID]})
hap_list=$(echo $(seq $base 16))

#chr=3
chr=$SLURM_ARRAY_TASK_ID
#hap_list=$(echo $(seq 2 16))
#pheno="male_flowering_d6"
#env="BLOIS_2014_OPT"

#for hapgrp in ${hap_list[@]}; do
#    Rscript GridLMM_pheno_x_env.R $pheno $env $chr $hapgrp 4
#    Rscript GridLMM_run.R $pheno $chr $hapgrp 4
#done


echo $SLURM_ARRAY_TASK_ID
while read pheno; do
    while read env; do
	for hapgrp in ${hap_list[@]}; do
	    echo $pheno
	    echo $env
	    echo $hapgrp
	    if [ $env == "ALL" ]
	    then
		#echo "ALL environments"
		Rscript GridLMM_run.R $pheno $chr $hapgrp 4
	    else
		#echo "Specific environment"
		Rscript GridLMM_pheno_x_env.R $pheno $env $chr $hapgrp 4
	    fi
	done
    done < env_list.txt
done < pheno_list.txt









