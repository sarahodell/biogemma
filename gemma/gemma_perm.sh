#!/bin/bash

rep=$1
pheno=$2
env=$3

gemma=/home/sodell/bin/./gemma-0.98.1-linux-static

cat phenotypes/${pheno}_${env}_phenotypes.txt | shuf > /scratch/sodell/${pheno}_${env}_phenotype_shuffle_rep${rep}.txt
phenofile=/scratch/sodell/${pheno}_${env}_phenotype_shuffle_rep${rep}.txt

for chr in {1..10}; do
  echo "Running $chr"
  genofile=/scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt
  #genofile=../genotypes/probabilities/allele_probs/test_bimbam.txt
  kinship=K_matrices/K_matrix_chr${chr}.txt
  outprefix=permute/GEMMA_wgs_GWAS_${pheno}_${env}_results_chr${chr}_perm_rep${rep}
  $gemma -g $genofile -p $phenofile -k $kinship -lmm -o $outprefix &
done

wait


echo "Finished rep $rep for $pheno $env successfully for all chromosomes" >> report.txt
#rm /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt
