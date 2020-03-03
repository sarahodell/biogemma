#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/tassel
#SBATCH -J pca
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem 16G

module load jdk/1.8
module load tassel

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx16g -fork1 -h Biogemma_600K_Genotypes_AGPv4.hmp.txt.gz -PrincipalComponentsPlugin -endPlugin -export PCA.txt

