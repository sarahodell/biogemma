#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J fillin
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 120:00:00

module load jdk/1.8

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx108g -FILLINFindHaplotypesPlugin -hmp hmpv3_founders/allchromosomes_hmp31_q30_founders.vcf.gz -o NAMHaps/ -extOut true -endPlugin
