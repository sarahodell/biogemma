#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J founders
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools/1.2
module load tabix

bgzip c6_hmp31_q30_edit.vcf
tabix -p c6_hmp31_q30_edit.vcf.gz
bcftools view -Oz -S nam_founders.txt c6_hmp31_q30_edit.vcf.gz > hmpv3_founders/c6_hmp31_q30_edit_founders.vcf.gz 


