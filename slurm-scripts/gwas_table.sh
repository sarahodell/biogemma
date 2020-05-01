#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM
#SBATCH -J gwas_table
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --array=1-7

module load R

declare -a env_list=("blank" "ALL" "BLOIS_2014_OPT" "BLOIS_2017_OPT" "GRANEROS_2015_OPT" "SZEGED_2017_OPT" "STPAUL_2017_WD" "NERAC_2016_WD")
env=$(echo ${env_list[$SLURM_ARRAY_TASK_ID]})

echo $env

#Plot haplotype probs
Rscript GridLMM_haplotypes/gwas_table.R $env

#Plot founder probs
Rscript  GridLMM_founderprobs/gwas_table.R $env

#Plot 600K 600K_SNP
Rscript  GridLMM_600KSNP/gwas_table.R $env
