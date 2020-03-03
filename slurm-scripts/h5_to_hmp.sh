#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J h5tohmp
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load jdk/1.8

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx16g -h5 ZeaGBSv27_publicSamples_rawGenos_AGPv3_20170206.h5 -FilterTaxaBuilderPlugin -taxaList NAM_RILS.json.gz -endPlugin -export NAMRILS_ZeaGBSv27_APGv3.hmp.txt.gz -exportType Hapmap

