#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/K_matrices
#SBATCH -J gridlmm
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=2
#SBATCH --mem=16G

module load R/3.5.2

#chr=10
chr=$SLURM_ARRAY_TASK_ID
hap_list=$(echo $(seq 2 16))

#echo $SLURM_ARRAY_TASK_ID
for hapgrp in $hap_list; do
    echo $hapgrp
    #Rscript make_ibd_kinship.R $chr $hapgrp
    Rscript ibd_K_founders.R $chr
done









