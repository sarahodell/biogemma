#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J founders
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load tabix/0.2.6
module load bcftools/1.2

for file in hmpv3_founders/*.vcf; do
    bgzip $file
    tabix -p vcf $file.gz
    bcftools view -Oz -S nam_founders.txt $file.gz > "${file%.*}"_founders.vcf.gz 
done


