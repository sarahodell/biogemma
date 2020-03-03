#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/uplifted_APGv4
#SBATCH -J hmpmerge
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --mem 32G

module load bcftools

bcftools merge --force-samples --file-list sample_hapmap.txt -Oz -o hmp321_agpv4_all.vcf.gz

    






