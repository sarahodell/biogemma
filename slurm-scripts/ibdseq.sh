#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/beagle
#SBATCH -J beagle
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --mem 30G
#SBATCH --ntasks=4

module load jdk


#java -Xmx32g -jar /home/sodell/bin/ibdseq.r1206.jar gt=Biogemma_WGS_Founders_filtered_snps_diploid_no_tester.vcf.gz out=Biogemma_WGS_Founders_IBDSEQ_2 r2max=0.95
java -Xmx29g -jar /home/sodell/bin/ibdseq.r1206.jar nthreads=4 chrom=$SLURM_ARRAY_TASK_ID gt=Biogemma_WGS_Founders_filtered_snps_diploid_no_tester.vcf.gz out=Biogemma_WGS_Founders_IBDSEQ_chr${SLURM_ARRAY_TASK_ID}
