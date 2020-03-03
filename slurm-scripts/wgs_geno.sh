#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/biogemma
#SBATCH -J fgeno
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem 32G
#SBATCH --ntasks=4

export PATH=/home/sodell/bin/bin:$PATH

#python make_diploid.py Biogemma_WGS_Founders_filtered_snps_reheader.vcf.gz Biogemma_WGS_Founders_filtered_snps_diploid.vcf

#python wgs_foundergeno2.py biogemma/Biogemma_Founders_600K_Genotypes_AGPv4_no_tester.vcf.gz Biogemma_0823_foundergeno

#python wgs_geno.py biogemma/WGS_positions.txt biogemma/Axiom600K_AGPv4.txt biogemma/Biogemma_DHLines_600K_Genotypes.vcf.gz Biogemma_WGS_geno

#python wgs_pmap.py biogemma/WGS_positions.txt biogemma/Axiom600K_AGPv4.txt Biogemma_WGS_pmap


#Rscript make_v4gmap.R

python wgs_founders.py Biogemma_WGS_Founders_filtered_snps_diploid_no_tester.vcf.gz Biogemma_WGS_alleles Biogemma_Founders_600K_Genotypes_AGPv4_no_tester.vcf.gz







