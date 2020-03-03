#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J hmpfilter
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load jdk/1.8
module load tabix
module load bcftools

gunzip c10_hmp321_withGL.vcf.gz
bgzip c10_hmp321_withGL.vcf
tabix -c vcf c10_hmp321_withGL.vcf.gz 
bcftools view -Oz -S nam_founders.txt c10_hmp321_withGL.vcf.gz > hmp321_founders_c10.vcf.gz 

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx16g -vcf hmp321_founders_c10.vcf.gz -HetsToUnknownPlugin -endPlugin -export hmp321_founders_noHets_c10.vcf.gz -exportType VCF

tabix -c vcf hmp321_founders_noHets_c10.vcf.gz

bcftools reheader -h c10_header.vcf -o hmp321_founders_c10_final.vcf.gz hmp321_founders_noHets_c10.vcf.gz

