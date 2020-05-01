#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J intersect
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --mem 7G

module load R
module load bedtools

Rscript scripts/match_to_bed.R $SLURM_ARRAY_TASK_ID

Rscript scripts/ibd_to_bed.R $SLURM_ARRAY_TASK_ID

Rscript scripts/ibd_to_bed2.R $SLURM_ARRAY_TASK_ID

Rscript scripts/bed_intersect.R $SLURM_ARRAY_TASK_ID

scripts/./bed_overlap.sh $SLURM_ARRAY_TASK_ID

Rscript scripts/overlap_quant.R $SLURM_ARRAY_TASK_ID
