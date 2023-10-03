#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim/qtl2_files
#SBATCH -J qtl2
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a_out.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a_error.txt
#SBATCH -t 24:00:00
#SBATCH --array=66-73
#SBATCH --ntasks=8
#SBATCH --mem=63G

module load R


#chr=${bob[]}
#block=${sue[$SLURM_ARRAY_TASK_ID]}

#json="Biogemma_WGS_c${chr}_${block}.json"
#echo $json
#Rscript qtl2_array.R $SLURM_ARRAY_TASK_ID

rep=$SLURM_ARRAY_TASK_ID
#chr=1
#chr=$SLURM_ARRAY_TASK_ID
for chr in {1..10}; do
  echo $chr
  Rscript ../qtl2_array.R $chr $rep 8
done
