#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J beagle
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=4
#SBATCH --mem 32G


module load jdk

chr=8
beagle=/home/sodell/bin/beagle.28Sep18.793.jar
data=biogemma/Biogemma_WGS_Founders_filtered_snps_diploid_no_tester.vcf.gz
donor=/scratch/sodell/hmp321_agpv4_chr8.vcf.gz
out=/scratch/sodell/Biogemma_WGS_Founders_chr${chr}_imputued.vcf.gz
map=beagle/plink_genetic.map

java -Xms20g -Xmx32g -jar $beagle gt=$data chrom=$chr impute=true  ap=true gp=true ref=$donor map=$map out=$out
