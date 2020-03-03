#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/MAGICSim_121718
#SBATCH -J buildmagic
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load bcftools

python ../impute/pscripts/build_simvcf.py MAGIC_DHSim_121318_breakpoints.txt MAGIC_DHSim_121718.vcf ../biogemma/founders --all True










