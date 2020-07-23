#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/gemma
#SBATCH -J gemma
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-500%50
#SBATCH --ntasks=8
#SBATCH --mem 63G

module load openblas
module load gsl
module load gemma

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" pheno_env_list_full.txt | cut -f3 -d,)"
#pheno="male_flowering_d6"
#env="BLOIS_2014_OPT"
#chr=$SLURM_ARRAY_TASK_ID
#chr=10
if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell ;
fi


if [ ! -f /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt ]; then
  cp ../genotypes/probabilities/allele_probs/bg${chr}_wgs_alleleprobs_bimbam.txt /scratch/sodell ;
fi
echo $pheno
echo $env
echo $chr

gemma=/home/sodell/bin/./gemma-0.98.1-linux-static
genofile=/scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt
#genofile=../genotypes/probabilities/allele_probs/test_bimbam.txt
phenofile=phenotypes/${pheno}_${env}_phenotypes.txt
kinship=K_matrices/K_matrix_chr${chr}.txt

$gemma -g $genofile -p $phenofile -k $kinship -lmm -o GEMMA_wgs_GWAS_${pheno}_${env}_results_chr${chr}



rm /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt
