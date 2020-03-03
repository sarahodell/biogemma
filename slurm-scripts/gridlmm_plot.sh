#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/GridLMM_600KSNP
#SBATCH -J glm_plot
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G

module load R/3.6.1

while read pheno; do
    echo $pheno
    while read env; do
	echo $env
        Rscript GridLMM_plot_by_chr.R $pheno $env
    done < env_list.txt
done < pheno_list.txt









