#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J rilmerge
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools
module load tabix
module load vcftools

for i in 1 2 3 4 5 6 7 8 9 10; do
    echo R${i} > R${i}.txt
    cat R${i}_B73xOh43_10samples.vcf | vcf-sort | bcftools reheader -s R${i}.txt -o R${i}_B73xOh43_10samples_edit.vcf
    bgzip R${i}_B73xOh43_10samples_edit.vcf
    tabix -p vcf R${i}_B73xOh43_10samples_edit.vcf.gz
    echo R${i}_B73xOh43_10samples_edit.vcf.gz >> ril_names.txt
done

bcftools merge -l ril_names.txt -m all -Oz > B73xOh43_RILSimAll_chr10.vcf.gz

    






