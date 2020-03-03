#!/bin/bash -l
module load R/3.6.0



for i in {1..500};do
    pheno="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f1 -d,)"
    env="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f2 -d,)"
    chr="$(sed "${i}q;d" pheno_env_list_full.txt | cut -f3 -d,)"

    echo $pheno
    echo $env
    echo $chr

    Rscript test.R $pheno $env $chr
done












