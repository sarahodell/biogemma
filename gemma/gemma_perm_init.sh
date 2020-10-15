#!/bin/bash

module load openblas
module load gsl
module load gemma

reps=3
pheno_envs=2
gemma=/home/sodell/bin/./gemma-0.98.1-linux-static

# Report file to check that all permutations run successfully
if [ -f report.txt ]; then
  rm report.txt
fi


#chr=$SLURM_ARRAY_TASK_ID
if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell ;
fi

for chr in {1..10}; do
  if [ ! -f /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt ]; then
    cp ../genotypes/probabilities/allele_probs/bg${chr}_wgs_alleleprobs_bimbam.txt /scratch/sodell ;
  fi
done


for pheno_env in $(seq 1 $pheno_envs); do
  pheno="$(sed "${pheno_env}q;d" pheno_env_list.txt | cut -f1 -d,)"
  env="$(sed "${pheno_env}q;d" pheno_env_list.txt | cut -f2 -d,)"
  echo $pheno
  echo $env
  for rep in $(seq 1 $reps); do
    #srun -t 1:00:00 ./test.sh &
    srun -t 4:00:00 ./gemma_perm.sh $rep $pheno $env &
    #pids[]{rep}]=$!
  done
done


#for pid in ${pids[*]}; do
#  wait $pid
#done
