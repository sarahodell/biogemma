#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J gbsfilter
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load jdk/1.8
module load tassel/5.2.14

for i in $(seq 1 10); do
    run_pipeline.pl -Xmx16g -h NAM_GBS/NAM_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr$i.hmp.txt.gz -HetsToUnknownPlugin -endPlugin -FilterTaxaPropertiesPlugin -minNotMissing 0.1 -endPlugin -export NAM_BPEC_AllZea_GBS_filtered_chr$i.hmp.txt.gz -exportType Hapmap
done


