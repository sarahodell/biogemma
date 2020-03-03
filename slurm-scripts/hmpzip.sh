#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/chr10hap
#SBATCH -J ziphmp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

for i in *; do
    bgzip -c $i > ../chr10hapzip/$i.gz 
done


