#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/MAGICSim_121318
#SBATCH -J builddh16
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools

python ../impute/pscripts/build_simvcf.py DH16Sim_121318.txt DH16_AGPv4_121318.vcf ../biogemma/founders --markerfile Axiom600K_chr10_AGPv4.bed










