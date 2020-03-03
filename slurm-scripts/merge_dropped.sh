#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J dropmerge
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools
module load tabix

for i in 1 2 3 4 5 6 7 8 9 10; do
    echo R${i} > RIL_sims/dropped/R${i}.txt
    awk -v OFS="\t" '$1=$1' RIL_sims/dropped/R${i}_B73xOh43_10samples_drop40_reheader.vcf | bcftools reheader -s R${i}.txt | bgzip -c > RIL_sims/dropped/R${i}_B73xOh43_10samples_drop40.vcf.gz
    tabix -p vcf RIL_sims/dropped/R${i}_B73xOh43_10samples_drop40.vcf.gz
    echo RIL_sims/dropped/R${i}_B73xOh43_10samples_drop40.vcf.gz >> RIL_sims/dropped/ril_names.txt
done

bcftools merge -l RIL_sims/dropped/ril_names.txt -m all -Oz > RIL_sims/dropped/B73xOh43_RILSimAll_chr10_drop40.vcf.gz

    






