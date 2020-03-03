#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/tmp
#SBATCH -J grabbam
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=4
#SBATCH --mem 32G

module load samtools

./bam_extract.sh

