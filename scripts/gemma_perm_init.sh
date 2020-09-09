#!/bin/bash

module load openblas
module load gsl
module load gemma

reps=1
pheno_envs=1
gemma=/home/sodell/bin/./gemma-0.98.1-linux-static



#chr=$SLURM_ARRAY_TASK_ID
if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell ;
fi

for chr in {1..10}; do
  if [ ! -f /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt ]; then
    cp genotypes/probabilities/allele_probs/bg${chr}_wgs_alleleprobs_bimbam.txt /scratch/sodell ;
  fi
done


for pheno_env in $(seq 1 $pheno_envs); do
  pheno="$(sed "${pheno_env}q;d" pheno_env_list.txt | cut -f1 -d,)"
  env="$(sed "${pheno_env}q;d" pheno_env_list.txt | cut -f2 -d,)"
  echo $pheno
  echo $env
  for rep in $(seq 1 $reps); do
    time scripts/gemma_perm.sh $rep $pheno $env &
    #pids[]{rep}]=$!
  done
done

wait

#for pid in ${pids[*]}; do
#  wait $pid
#done
