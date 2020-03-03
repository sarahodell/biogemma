#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J vcf_edit
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

awk '/^#CHROM/ {printf("##contig=<ID=1,length=301476924>\n");} {print;}' /home/sodell/projects/impute/hmpv3/c1_hmp31_q30.vcf > c1_hmp31_q30_edit.vcf


