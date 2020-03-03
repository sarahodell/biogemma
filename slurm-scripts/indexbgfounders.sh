#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/biogemma
#SBATCH -J founders
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load tabix/0.2.6
module load bcftools/1.2

#bcftools reheader -h founderheader_v4.vcf -o Biogemma_600K_Genotypes_AGPv4.vcf Biogemma_600K_Genotypes_AGPv4.vcf
#bgzip Biogemma_600K_Genotypes_AGPv4.vcf
#tabix -p vcf Biogemma_600K_Genotypes.vcf.gz

#bcftools view -Oz -S founders/biogemma_founders.txt Biogemma_600K_Genotypes_AGPv4.vcf.gz > BiogemmaFounders_600K_Genotypes_AGPv4.vcf.gz
tabix -p vcf BiogemmaFounders_600K_Genotypes_AGPv4.vcf.gz

while read line; do
    sample=$line
    bcftools view -s $sample -Ov BiogemmaFounders_600K_Genotypes_AGPv4.vcf.gz > founders/${sample}_600K_Genotypes_AGPv4.vcf
    bgzip founders/${sample}_600K_Genotypes_AGPv4.vcf
    tabix -p vcf founders/${sample}_600K_Genotypes_AGPv4.vcf.gz
done < founders/biogemma_founders.txt
