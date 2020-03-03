#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J gbsfilter
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load jdk/1.8

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx16g -h NAMRILS_ZeaGBSv27_APGv3.hmp.txt -HetsToUnknownPlugin -endPlugin -export NAMRILS_ZeaGBSv27_APGv3_noHets2.hmp.txt.gz -exportType Hapmap


