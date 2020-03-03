#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J vcf_edit10
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00


gunzip -k /home/sodell/projects/impute/c10_hmp321_withGL.vcf.gz
awk '/^#CHROM/ {printf("##contig=<ID=10,length=149632204>\n");} {print;}' /home/sodell/projects/impute/c10_hmp321_withGL.vcf > c10_hmp321_reheader.vcf


