#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/
#SBATCH -J fillin
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --mem=125G
#SBATCH --ntasks=16

module load jdk
module load tassel

#for i in $(seq 1 10); do
#   run_pipeline.pl -Xmx32g -h NA#M_GBS/NAM_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr$i.hmp.txt.gz -HetsToUnknownPlugin -endPlugin -FilterTaxaPropertiesPlugin -minNotMissing 0.1 -endPlugin -export NAM_BPEC_AllZea_GBS_filtered_chr$i.hmp.txt.gz -exportType Hapmap
#done


run_pipeline.pl -Xms64g -Xmx84g -FILLINFindHaplotypesPlugin -hmp /scratch/sodell/hmp321_agpv4_chr8.vcf.gz -nV true -o /scratch/sodell/HapMap3_chr8

#run_pipeline.pl -Xmx16g -FILLINImputationPlugin -hmp FILLIN/Biogemma_WGS_Founders_chr8.vcf.gz -d HapMap3_chr8 -o Biogemma_WGS_Founders_chr8_imputed

#run_pipeline.pl -Xms80g -Xmx120g -vcf /scratch/sodell/hmp321_agpv4_chr8.vcf.gz -ImputationPlugin -nearestNeighbors 5 -endPlugin -export /scratch/sodell/hmp321_agpv4_chr8_imputed.vcf.gz -exportType VCF



