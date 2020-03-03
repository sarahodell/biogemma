#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J crossmap
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00


python ./CrossMap-0.2.5/bin/CrossMap.py vcf AGPv2_to_AGPv3.chain.gz maize600k_eliteLines.vcf ../AGPv3/Zea_mays.AGPv3.22.dna.genome.fa maize600k_eliteLines_v3.vcf












