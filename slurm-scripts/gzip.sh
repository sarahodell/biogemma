#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/WGS
#SBATCH -J zipzip
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --array=1-10
#SBATCH --mem 3G

gzip Biogemma_WGS_all_alleles_final_chr${SLURM_ARRAY_TASK_ID}.txt

