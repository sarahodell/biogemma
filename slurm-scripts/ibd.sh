#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J haploprobs
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=8
#SBATCH --mem 60G

module load R

chr=$SLURM_ARRAY_TASK_ID
#Rscript scripts/germline_haplotypes.R $chr
#if [ -f genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt ]; then
#  rm genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt
#fi
#base_list=(blank 8 7 7 8 6 7 7 8 7 7)
#base=$(echo ${base_list[$chr]})
#hap_list=$(echo $(seq $base 16))

#Rscript scripts/build_ibd_segments.R $chr
base_list=( blank $(cut -d$'\t' -f2 genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt) )
base=$(echo ${base_list[$chr]})
hap_list=$(echo $(seq $base 16))

echo "Built IBD Segments"
for hapgrp in ${hap_list[@]}; do
  echo "starting haplogroup $hapgrp"
  Rscript scripts/build_haplotypes.R $chr $hapgrp
  echo "Created Haplotype probabilities"
  Rscript scripts/filter_haplotypes.R $chr $hapgrp
  echo "Filtered Haplotype probabilities"
  #Rscript scripts/filter_haplotypes.R $chr $hapgrp &
  #echo "Filtered Haplotype probabilities"
done
