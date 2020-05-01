#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/vgt1
#SBATCH -J samtools
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem=15G
#SBATCH --ntasks=2

module load samtools

scripts/./bam_extract.sh












