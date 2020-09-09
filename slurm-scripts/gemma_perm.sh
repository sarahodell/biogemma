#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J gemma_perm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -w bigmem3
#SBATCH -t 96:00:00
#SBATCH --ntasks=10
#SBATCH --mem 78G

module load openblas
module load gsl
module load gemma

rep=1
gemma=/home/sodell/bin/./gemma-0.98.1-linux-static
#pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" rerun_pheno_env_list.txt | cut -f1 -d,)"
#env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" rerun_pheno_env_list.txt | cut -f2 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" rerun_pheno_env_list.txt | cut -f3 -d,)"
pheno="grain_yield_15"
env="NERAC_2016_WD"

#chr=$SLURM_ARRAY_TASK_ID
if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell ;
fi

cat gemma/phenotypes/${pheno}_${env}_phenotypes.txt | shuf > /scratch/sodell/${pheno}_${env}_phenotype_shuffle.txt
phenofile=/scratch/sodell/${pheno}_${env}_phenotype_shuffle.txt
#echo $pheno
#echo $env
#echo $chr
for chr in {1..10}; do
  if [ ! -f /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt ]; then
    cp genotypes/probabilities/allele_probs/bg${chr}_wgs_alleleprobs_bimbam.txt /scratch/sodell ;
  fi
  genofile=/scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt
  #genofile=../genotypes/probabilities/allele_probs/test_bimbam.txt
  kinship=gemma/K_matrices/K_matrix_chr${chr}.txt
  outprefix=/scratch/sodell/GEMMA_wgs_GWAS_${pheno}_${env}_results_chr${chr}_perm_rep${rep}
  srun -n 1 $($gemma -g $genofile -p $phenofile -k $kinship -lmm -o $outprefix) &
  echo "Running $chr"
done

wait

#rm /scratch/sodell/bg${chr}_wgs_alleleprobs_bimbam.txt
