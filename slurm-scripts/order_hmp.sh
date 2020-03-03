#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J orderhmp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load jdk/1.8
module load tassel/5.2.14


run_pipeline.pl -Xmx16G -SortGenotypeFilePlugin -inputFile biogemma/biogemma_600K_Genotypes_AGPv4.hmp.txt -outputFile biogemma/Biogemma_600K_Genotypes_AGPv4.hmp.txt -fileType Hapmap



