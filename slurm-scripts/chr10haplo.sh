#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J hap10
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 120:00:00

module load jdk/1.8

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx32g -FILLINFindHaplotypesPlugin -hmp hmp3_founders2/hmp3_founders_final.vcf.gz -o chr10hap2 -minTaxa 1 -mxErr 0 -extOut true -endPlugin
