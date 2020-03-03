#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/
#SBATCH -J founderfile
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH --array=1-10%5
#SBATCH -t 24:00:00

module load bcftools

bcftools view -Oz -r $SLURM_ARRAY_TASK_ID biogemma/Biogemma_DHLines_600K_Genotypes.vcf.gz > biogemma/Biogemma_DHLines_600K_Genotypes_chr$SLURM_ARRAY_TASK_ID.vcf.gz

python impute/qtl2/genofile2.py biogemma/Biogemma_DHLines_600K_Genotypes_chr$SLURM_ARRAY_TASK_ID.vcf.gz Biogemma_DHgenos/DH_geno_chr${SLURM_ARRAY_TASK_ID}_121718.csv












