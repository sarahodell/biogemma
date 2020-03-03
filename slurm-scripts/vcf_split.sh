#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J split
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools/1.2

for sample in $(bcftools query -l hmpv3_founders/allchromosomes_hmp31_q30_founders.vcf.gz); do
    echo $sample > hmpv3_founders/$sample.txt
    bcftools view -Oz -S hmpv3_founders/$sample.txt hmpv3_founders/allchromosomes_hmp31_q30_founders.vcf.gz > hmpv3_founders/${sample}_hmp321_all.vcf.gz
done
    






