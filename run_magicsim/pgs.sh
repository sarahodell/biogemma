#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J pgs
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-101
#SBATCH --mem 11G
#SBATCH --ntasks 3

#Creates a merged vcf files of 400 Simulated MAGIC lines

module load R

rep=$SLURM_ARRAY_TASK_ID

#Rscript rrblup_sim.R $rep
Rscript alt_pgs.R $rep
