#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J beagle
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=16
#SBATCH --mem 125G


module load jdk

chr=1
beagle=/home/sodell/bin/beagle.28Sep18.793.jar
data=/scratch/sodell/Biogemma_WGS_Founders_filtered_snps_diploid_phased.vcf.gz
#donor=/scratch/sodell/hmp321_agpv4_chr${chr}.vcf.gz
out=/scratch/sodell/Biogemma_WGS_Founders_chr${chr}_imputed
map=genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr${chr}.map

java -Xmx125g -jar $beagle gt=$data chrom=$chr nthreads=16 impute=true ap=true gp=true map=$map out=$out
