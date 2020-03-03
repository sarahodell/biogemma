#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_050719
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=8
#SBATCH --mem=64G

module load R/3.5.2


#chr=${bob[$SLURM_ARRAY_TASK_ID]}
#block=${sue[$SLURM_ARRAY_TASK_ID]}

#json="Biogemma_WGS_c${chr}_${block}.json"
#echo $json
#Rscript qtl2_array.R $SLURM_ARRAY_TASK_ID

echo $SLURM_ARRAY_TASK_ID
#Rscript Biogemma_121318/wgs_alleleprobs.R $SLURM_ARRAY_TASK_ID
Rscript qtl2_array.R $SLURM_ARRAY_TASK_ID










