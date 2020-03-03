#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=50
#SBATCH --mem=128G

module load R/3.3.1

#bob=( 1 1 1 1 )
#sue=( 1 2 3 4 )

#chr=${bob[$SLURM_ARRAY_TASK_ID]}
#block=${sue[$SLURM_ARRAY_TASK_ID]}

#json="Biogemma_WGS_c${chr}_${block}.json"
#echo $json
#Rscript qtl2_array.R 1 1
Rscript qtl2_demo.R











