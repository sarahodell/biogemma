#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J genofiles
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-100%10
#SBATCH --mem 7G
#SBATCH --ntasks 1

#Creates a qtl2 geno file of 400 Simulated MAGIC lines
module load python/2.7.15
module load bcftools

rep=$SLURM_ARRAY_TASK_ID
#rep=74
#25 39 74
tabix -p vcf merged_vcfs/MAGIC_DHSimAll_rep${rep}_v2.vcf.gz
for c in {1..10}; do
  python genofile2.py merged_vcfs/MAGIC_DHSimAll_rep${rep}_v2.vcf.gz qtl2_files/MAGIC_DHSimAll_rep${rep}_chr${c}_v2.csv --region $c
done
