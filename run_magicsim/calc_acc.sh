#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J calc_acc
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 96:00:00
#SBATCH --array=1-100
#SBATCH --ntasks=4
#SBATCH --mem=15G

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
  Rscript compare_props.R $rep $chr 4
done
