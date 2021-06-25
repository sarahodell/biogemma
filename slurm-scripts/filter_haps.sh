#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J corfilter
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a_out-.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a_error.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-100
#SBATCH --ntasks=2
#SBATCH --mem=7G

module load R

rep=$SLURM_ARRAY_TASK_ID
#chr=10

#echo $SLURM_ARRAY_TASK_ID
#Rscript breakup_haplos.R $SLURM_ARRAY_TASK_ID
#Rscript scripts/filter_geno.R $chr
#Rscript breakup_haplos.R $chr
for c in {1..10};do
  Rscript filter_geno.R $c $rep
done
