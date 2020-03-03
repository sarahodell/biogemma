#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM
#SBATCH -J gridlmm
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load conda3
module load R/3.5.2
source activate myRenv3_5

pheno="male_flowering_days"
#env="ALL"
env="BLOIS_2014_OPT"
chr=$SLURM_ARRAY_TASK_ID
hap_list=$(echo $(seq 2 16))

#echo $SLURM_ARRAY_TASK_ID
for hapgrp in $hap_list; do
    echo $pheno
    echo $env
    echo $hapgrp
    Rscript GridLMM_run.R $pheno $chr $hapgrp 4
    #Rscript GridLMM_pheno_x_env.R $pheno $env $chr $hapgrp 4
done

Rscript GridLMM_plot.R $pheno $env $chr 4








