#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/GridLMM_600KSNP
#SBATCH -J glm_snp
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=4
#SBATCH --mem=4G
#SBATCH --array=1-50%10

module load R

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f2 -d,)"

echo $pheno
echo $env

for c in {1..10};do
    if [ $env == "ALL" ]
    then
	echo "ALL environments"
	Rscript GridLMM_600K_blup_QTLbounds.R $pheno $c 4
    else
	echo "Specific environment"
	Rscript GridLMM_600K_QTLbounds.R $pheno $env $c 4
    fi
done








