#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/MAGICSim_121718
#SBATCH -J rilmerge
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools
module load tabix
module load vcftools

for i in {1..400}; do
    echo Sim${i} > Sim${i}.txt
    cat Sim${i}_MAGIC_DHSim_121718.vcf | vcf-sort | bcftools reheader -s Sim${i}.txt -o Sim${i}_MAGIC_DHSim_121718_edit.vcf
    bgzip Sim${i}_MAGIC_DHSim_121718_edit.vcf
    tabix -p vcf Sim${i}_MAGIC_DHSim_121718_edit.vcf.gz
    echo Sim${i}_MAGIC_DHSim_121718_edit.vcf.gz >> line_names.txt
done

bcftools merge -l line_names.txt -m all -Oz > MAGIC_DHSimAll_121718_chr10.vcf.gz

    






