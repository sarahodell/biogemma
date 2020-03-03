#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/hmpv3_founders
#SBATCH -J concat
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools

bcftools concat -Oz c1_hmp31_q30_edit_founders.vcf.gz c2_hmp31_q30_edit_founders.vcf.gz c3_hmp31_q30_edit_founders.vcf.gz c4_hmp31_q30_edit_founders.vcf.gz c5_hmp31_q30_edit_founders.vcf.gz c6_hmp31_q30_edit_founders.vcf.gz c7_hmp31_q30_edit_founders.vcf.gz c8_hmp31_q30_edit_founders.vcf.gz c9_hmp31_q30_edit_founders.vcf.gz c10_hmp31_q30_edit_founders.vcf.gz > allchromosomes_hmp31_q30_founders.vcf.gz
    






