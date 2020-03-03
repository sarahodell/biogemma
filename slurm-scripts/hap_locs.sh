#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/chr10hap2
#SBATCH -J haplocs
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

echo -e "haplotype\tstart\tend" > hap_locs.txt
for file in *; do
    name=$(echo $file | cut -d'.' -f2)
    start=$(sed -n '2{p;q;}' $file | cut -f4) 
    end=$(tail -n 1 $file | cut -f4)
    echo -e "$name\t$start\t$end" >> hap_locs.txt
done
