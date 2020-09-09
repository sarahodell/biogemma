#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM
#SBATCH -J glm_plot
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --array=1-50%10

module load R

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list.txt | cut -f2 -d,)"

echo $pheno
echo $env

#Plot haplotype probs
Rscript GridLMM_haplotypes/GridLMM_plot_by_chr.R $pheno $env

#Plot founder probs
Rscript  GridLMM_founderprobs/GridLMM_plot_by_chr.R $pheno $env

#Plot 600K 600K_SNP
Rscript  GridLMM_600KSNP/GridLMM_plot_by_chr.R $pheno $env

echo "Finished array task $SLURM_ARRAY_TASK_ID" >> report.txt
