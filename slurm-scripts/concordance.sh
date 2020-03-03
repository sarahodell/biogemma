#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/biogemma
#SBATCH -J conc
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa
bcf=/home/sodell/bin/bin/bcftools

#$bcf convert -Oz --haploid2diploid Biogemma_WGS_Founders_600K_Genotypes.vcf.gz > Biogemma_WGS_Founders_600K_Genotypes_diploid.vcf.gz
tabix -p vcf Biogemma_WGS_Founders_600K_Genotypes_diploid.vcf.gz

java -Xmx16g -jar $GATK \
     -T GenotypeConcordance \
     -R $genome \
     -eval Biogemma_WGS_Founders_600K_Genotypes_diploid.vcf.gz \
     -comp Biogemma_Founders_600K_Genotypes_AGPv4.vcf.gz \
     -o WGS_vs_600K_output.grp












