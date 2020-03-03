#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/plink_files
#SBATCH -J haploprobs
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=3
#SBATCH --mem 23G

module load R

#Rscript qtl2.R 
#Rscript sim_geno.R
#Rscript bg_percentage.R
#Rscript control_file.R
#Rscript get_ibd.R $SLURM_ARRAY_TASK_ID
#Rscript ibd_segments.R $SLURM_ARRAY_TASK_ID
#Rscript wgs_genoprobs.R
#Rscript recode_ibd.R
#Rscript bg_coverage.R 
#Rscript ibd_blocks.R $SLURM_ARRAY_TASK_ID
#Rscript bg_coverage.R $SLURM_ARRAY_TASK_ID

#IBDFILE="ibd_segments/bg${SLURM_ARRAY_TASK_ID}_wgs_ibd_segments_adjusted.txt"
#GENOFILE="geno_probs/bg${SLURM_ARRAY_TASK_ID}_genoprobs_010319.rds"
#PMAP="qtl2_startfiles/Biogemma_pmap_c${SLURM_ARRAY_TASK_ID}.csv"
#OUTFILE="haplotype_probs/bg${SLURM_ARRAY_TASK_ID}_haplotype_probs_012720.txt"

#Rscript Rscripts/ibd_blocks.R $SLURM_ARRAY_TASK_ID $IBDFILE $GENOFILE $PMAP $OUTFILE

Rscript germline_haplotypes.R 10
