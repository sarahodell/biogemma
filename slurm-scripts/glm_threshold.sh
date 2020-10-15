#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM
#SBATCH -J glm_threshold
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R

if [ -f threshold_table.txt ]; then
  rm threshold_table.txt
fi


for i in {1..50}; do
# Permuation 1000 times
  pheno="$(sed "${i}q;d" pheno_env_list.txt | cut -f1 -d,)"
  env="$(sed "${i}q;d" pheno_env_list.txt | cut -f2 -d,)"
  #pheno="harvest_grain_moisture"
  #env="BLOIS_2014_OPT"
  echo $pheno
  echo $env

  #Haplotype threshold
  Rscript GridLMM_haplotypes/permute_haplo/find_threshold.R $pheno $env
  wait
  #Founder threshold
  Rscript GridLMM_founderprobs/permute/find_threshold.R $pheno $env
  wait
  #600K threshold
  Rscript GridLMM_600KSNP/permute/find_threshold.R $pheno $env
  wait
done

echo "Finished array task $SLURM_ARRAY_TASK_ID" >> report.txt
