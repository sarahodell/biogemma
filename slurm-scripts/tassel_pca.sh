#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J pca
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem 16G

module load jdk/1.8
module load tassel

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx8g -vcf genotypes/600K/Biogemma_Founders_600K_Genotypes_AGPv4_sorted.vcf.gz -tree Neighbor -treeSaveDistance false -export genotypes/600K/tree.nj.txt
