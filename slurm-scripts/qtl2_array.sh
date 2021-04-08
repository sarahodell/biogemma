#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a_out.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a_error.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=8
#SBATCH --mem=65G

module load R


#chr=${bob[$SLURM_ARRAY_TASK_ID]}
#block=${sue[$SLURM_ARRAY_TASK_ID]}

#json="Biogemma_WGS_c${chr}_${block}.json"
#echo $json
#Rscript qtl2_array.R $SLURM_ARRAY_TASK_ID
chr=1
#chr=$SLURM_ARRAY_TASK_ID
echo $chr
#Rscript Biogemma_121318/wgs_alleleprobs.R $SLURM_ARRAY_TASK_ID
Rscript scripts/qtl2_array.R $chr 8
