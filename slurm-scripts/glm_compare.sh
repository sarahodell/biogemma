#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM
#SBATCH -J glm_compare
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R/3.6.0


chr=$SLURM_ARRAY_TASK_ID
#hap_list=$(echo $(seq 2 16))

pheno="male_flowering_d6"

while read env; do
    for hapgrp in $hap_list; do
	echo $env
        echo $hapgrp
    done
done < env_list.txt










