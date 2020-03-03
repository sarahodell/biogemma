#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J founders
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load tabix/0.2.6
module load bcftools/1.2

while read line; do
    sample=$line
    tabix -p vcf hmp3_founders2/${sample}_c10_hmp321_final.vcf.gz    
done < nam_founders.txt


