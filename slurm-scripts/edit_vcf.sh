#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J vcf_edit2
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00


gunzip -k /home/sodell/projects/impute/hmpv3/c2_hmp31_q30.vcf.gz
awk '/^#CHROM/ {printf("##contig=<ID=2,length=237917468>\n");} {print;}' /home/sodell/projects/impute/hmpv3/c2_hmp31_q30.vcf > c2_hmp31_q30_edit.vcf


