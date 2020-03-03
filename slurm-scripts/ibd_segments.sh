#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J all
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --array=1-10%3
#SBATCH --ntasks=5
#SBATCH --mem=40G

module load R/3.5.1

#bob=( 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 )
#sue=( 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 )

#chr=${bob[$SLURM_ARRAY_TASK_ID]}
#block=${sue[$SLURM_ARRAY_TASK_ID]}

#json="Biogemma_WGS_c${chr}_${block}.json"
#echo $json
#Rscript qtl2_array.R $chr $block

echo $SLURM_ARRAY_TASK_ID
Rscript wgs_founders3.R $SLURM_ARRAY_TASK_ID











