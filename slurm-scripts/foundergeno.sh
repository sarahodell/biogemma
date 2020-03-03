#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/
#SBATCH -J founderfile
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH --array=1-10%5
#SBATCH -t 24:00:00

module load bcftools

bcftools view -Oz -r $SLURM_ARRAY_TASK_ID biogemma/Biogemma_Founders_600K_Genotypes_AGPv4_no_tester.vcf.gz > biogemma/Biogemma_Founders_600K_Genotypes_AGPv4_no_tester_chr${SLURM_ARRAY_TASK_ID}.vcf.gz

python impute/qtl2/foundergeno.py biogemma/Biogemma_Founders_600K_Genotypes_AGPv4_no_tester_chr${SLURM_ARRAY_TASK_ID}.vcf.gz Biogemma_foundergenos/Founder_genos_chr${SLURM_ARRAY_TASK_ID}_121718.csv












