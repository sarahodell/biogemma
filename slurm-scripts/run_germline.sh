#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J germline
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem=15G
#SBATCH --ntasks=4

GERMLINE="/home/sodell/bin/germline-1-5-3/./germline"
PEDFILE="plink_files/Biogemma_Founders_600K_chr10.ped"
MAPFILE="plink_files/Biogemma_Founders_600K_chr10.map"
#GMAPFILE="plink_files/Biogemma_Founders_600K_chr10_gmap.map"
OUTPUT="plinkfiles/Biogemma_germline_IBD_chr10"

$GERMLINE -map $GMAPFILE -input $PEDFILE $MAPFILE -output $OUTPUT



