#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J h5schema
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

run_pipeline.pl -hdf5Schema ZeaGBSv27_publicSamples_rawGenos_AGPv3_20170206.h5 -export rawGenos_h5_schema.txt

