#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J gbs10
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools/1.2


bcftools view -Oz -r 10 NAMRILS_ZeaGBSv27_APGv3_noHets.hmp.txt.gz > NAMRILS_ZeaGBSv27_AGPv3_chr10

    






