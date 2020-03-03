#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/Biogemma_121318/GridLMM/GridLMM_founderprobs/no_filtering
#SBATCH -J glm_nofilt
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-10%5
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R/3.5.2

#chr=3
chr=$SLURM_ARRAY_TASK_ID
hap_list=$(echo $(seq 2 16))

#echo $SLURM_ARRAY_TASK_ID
while read pheno; do
    while read env; do
	for hapgrp in $hap_list; do
	    echo $pheno
	    echo $env
	    echo $hapgrp
	    if [ $env == "ALL" ]
	    then
		#echo "ALL environments"
		Rscript GridLMM_run.R $pheno $chr $hapgrp 4
	    else
		#echo "Specific environment"
		Rscript GridLMM_pheno_x_env.R $pheno $env $chr $hapgrp 4
	    fi
	done
	echo "Building Plots"
#        Rscript GridLMM_plot.R $pheno $env $chr 4
#        Rscript GridLMM_sig_filter.R $pheno $env $chr
    done < env_list.txt
done < pheno_list.txt









