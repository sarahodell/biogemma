#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/GridLMM_founderprobs
#SBATCH -J glm_fp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=3G

module load R

#pheno="asi"

chr=$SLURM_ARRAY_TASK_ID

while read pheno; do
    echo $pheno
    while read env; do
        if [ $env == "ALL" ]
        then
            echo "ALL environments"
            Rscript GridLMM_run_founders.R $pheno $chr 4
        else
            echo $env
            Rscript GridLMM_pheno_x_env.R $pheno $env $chr 4
        fi
    done < env_list.txt
done < pheno_list.txt










