#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J genofile
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
module load bcftools

#header='id,chr,pos,A632_usa,B73_inra,CO255_inra,FV252_inra,OH43_inra,A654_inra,FV2_inra,C103_inra,EP1_inra,D105_inra,W117_inra,B96,DK63,F492,ND245,VA85'

#echo $header > biogemma/BiogemmaFounders.csv
#bcftools query -f '%ID,%CHROM,%POS[,%GT]\n' biogemma/BiogemmaFounders_600K_Genotypes.vcf.gz >> biogemma/BiogemmaFounders.csv

python pscripts/foundergeno.py biogemma/BiogemmaFounders.csv qtl2/MAGICSim_062018/MAGICSim_foundergeno.csv










