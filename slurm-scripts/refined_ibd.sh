#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J refined_ibd
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=1
#SBATCH --mem 8G


module load jdk

chr=$SLURM_ARRAY_TASK_ID
refinedibd=/home/sodell/bin/refined-ibd.17Jan20.102.jar
data=biogemma/Biogemma_WGS_Founders_no_missing_phased.vcf.gz
out=beagle/Biogemma_WGS_Founders_IBD_chr${chr}
map=beagle/plink_genetic.map

java -Xss5m -Xmx7g -jar $refinedibd nthreads=1 gt=$data window=10.0 chrom=$chr map=$map out=$out
