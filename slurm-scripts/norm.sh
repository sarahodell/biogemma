#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J norm
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

export PATH=/home/sodell/bin/bin/:$PATH

bcftools norm -Oz --do-not-normalize --check-ref ws --fasta-ref /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa --output biogemma/Biogemma_600K_Genotypes_AGPv4_norm_final.vcf biogemma/Biogemma_600K_Genotypes_AGPv4_final.vcf.gz











