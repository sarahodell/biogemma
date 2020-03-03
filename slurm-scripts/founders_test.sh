#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J founders
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load tabix/0.2.6
module load bcftools/1.2

bgzip c2_hmp31_q30_edit.vcf
tabix -p vcf c2_hmp31_q30_edit.vcf.gz
bcftools view -Oz -S nam_founders.txt c2_hmp31_q30_edit.vcf.gz > c2_hmp31_q30_edit_founders.vcf.gz
