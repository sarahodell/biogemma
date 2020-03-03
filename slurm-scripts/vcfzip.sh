#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J zipvcf
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load tabix/0.2.6

for chr in 1 3 4 5 6 7 8 9 10; do
    VCF=hmpv3_founders/c${chr}_hmp31_q30_edit.vcf
    echo $VCF
    bgzip $VCF
    tabix -p vcf $VCF.gz   
done


