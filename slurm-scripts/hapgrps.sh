#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/ibd_121018
#SBATCH -J hapgrp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=32G

module load R/3.5.1


echo $SLURM_ARRAY_TASK_ID
Rscript haplogroups.R $SLURM_ARRAY_TASK_ID











