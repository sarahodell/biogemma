#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/MITE_only
#SBATCH -J glm_fp
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=3G

module load R

#pheno="asi"
pheno="male_flowering_d6"
chr=$SLURM_ARRAY_TASK_ID


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
