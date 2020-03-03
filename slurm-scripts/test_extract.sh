#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J founders
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools/1.2

for i in 1 3 4 5 6 7 8 9 10; do
    bcftools view -Oz -S nam_founders.txt hmpv3_founders/c${i}_hmp31_q30_edit.vcf.gz > hmpv3_founders/c${i}_hmp31_q30_edit_founders.vcf.gz 
done

