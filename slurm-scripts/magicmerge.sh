#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/MAGICSim_062018
#SBATCH -J magicmerge
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools
module load tabix
module load vcftools

for i in 1 2 3 4 5 6 7 8 9 10; do
    echo M${i} > M${i}.txt
    cat M${i}_MAGICSim_Biogemma_062018.vcf | vcf-sort | bcftools reheader -s M${i}.txt -o M${i}_MAGICSim_Biogemma_062018_edit.vcf
    bgzip M${i}_MAGICSim_Biogemma_062018_edit.vcf
    tabix -p vcf M${i}_MAGICSim_Biogemma_062018_edit.vcf.gz
    echo M${i}_MAGICSim_Biogemma_062018_edit.vcf.gz >> line_names.txt
done

bcftools merge -l line_names.txt -m all -Oz > MAGICSim_Biogemma_062018.vcf.gz
bcftools reheader -h ../MAGICheader.vcf -o MAGICSim_BiogemmaFull_062018.vcf.gz MAGICSim_Biogemma_062018.vcf.gz
    






