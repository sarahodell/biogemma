#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J refined_ibd
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem 31G


module load jdk

chr=$SLURM_ARRAY_TASK_ID
refinedibd=/home/sodell/bin/refined-ibd.17Jan20.102.jar
data=genotypes/plink_files/WGS/Biogemma_WGS_Founders_no_missing_biallelic.vcf.gz
out=ibd_segments/refinedibd/Biogemma_WGS_Founders_RefinedIBD_chr${chr}
map=ibd_segments/plink_genetic.map

java -Xss5m -Xmx30g -jar $refinedibd nthreads=4 gt=$data window=10.0 length=0.2 chrom=$chr map=$map out=$out
