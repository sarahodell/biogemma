#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J alleles
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=6
#SBATCH --mem=48G

module load R

#echo $SLURM_ARRAY_TASK_ID
#chr=$SLURM_ARRAY_TASK_ID
chr=2
#Rscript scripts/wgs_alleleprobs.R $SLURM_ARRAY_TASK_ID
if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell ;
fi


if [ ! -f /scratch/sodell/Biogemma_WGS_all_alleles_final_chr${chr}.txt.gz ]; then
  cp genotypes/WGS/Biogemma_WGS_all_alleles_final_chr${chr}.txt.gz /scratch/sodell ;
fi

Rscript scripts/wgs_alleleprobs.R $chr

rm /scratch/sodell/Biogemma_WGS_all_alleles_final_chr${chr}.txt.gz
