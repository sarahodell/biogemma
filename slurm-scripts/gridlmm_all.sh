#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/GridLMM
#SBATCH -J gridlmm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=4
#SBATCH --mem=14G

module load R
cd GridLMM_600KSNP
echo "600K_SNP"
for i in {1..570}; do
  pheno="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
  env="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
  chr="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f3 -d,)"
  echo $pheno
  echo $chr
  if [ $env == "ALL" ]
  then
    echo "ALL environments"
    Rscript GridLMM_600K_all.R $pheno $chr 4
  else
    echo $env
    Rscript GridLMM_600K_pheno_x_env.R $pheno $env $chr 4
  fi
done

cd ..
cd GridLMM_founderprobs
echo "Founder_probs"
for i in {1..570}; do
  pheno="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
  env="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
  chr="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f3 -d,)"

  echo $pheno
  echo $chr

  if [ $env == "ALL" ]
  then
    echo "ALL environments"
    Rscript GridLMM_run_founders.R $pheno $chr 4
  else
    echo $env
    Rscript GridLMM_pheno_x_env.R $pheno $env $chr 4
  fi
done

cd ..
cd GridLMM_haplotypes
echo "Haplotype_probs"
for i in {1..570}; do
  pheno="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
  env="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
  chr="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f3 -d,)"

  base_list=( blank $(cut -d$'\t' -f2 ../../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt) )
  base=$(echo ${base_list[$chr]})
  hap_list=$(echo $(seq $base 16))
  echo $chr
  echo $hap_list
  #base=$(echo ${base_list[$chr]})
  #hap_list=$(echo $(seq $base 16))

  for hapgrp in ${hap_list[@]}; do
      echo $pheno
      echo $env
      echo $chr
      echo $hapgrp
      if [ $env == "ALL" ]
      then
  	     Rscript GridLMM_run.R $pheno $chr $hapgrp 4
      else
  	     Rscript GridLMM_pheno_x_env.R $pheno $env $chr $hapgrp 4
      fi
      #echo "$pheno $env $chr $hapgrp done" >> report.txt
  done
done
