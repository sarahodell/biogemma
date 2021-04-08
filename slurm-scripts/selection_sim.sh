#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J chisq
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-100
#SBATCH --ntasks=3
#SBATCH --mem 11G

module load R

rep=$SLURM_ARRAY_TASK_ID

for i in {1..10}; do
  #Rscript filter_geno.R $i $rep
  Rscript founder_rep.R $i $rep
done
