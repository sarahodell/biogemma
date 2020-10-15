#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM/MITE_only
#SBATCH -J glm_plot
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G

module load R

pheno="male_flowering_d6"


while read env; do
    Rscript GridLMM_plot_by_chr.R $pheno $env
done < env_list.txt
