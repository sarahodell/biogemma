#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J markers
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
module load bcftools

python marker_generator.py hmp3_founders2/hmp3_founders_final.vcf.gz c10_markers.txt 500000 










