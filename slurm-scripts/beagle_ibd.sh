#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J beagle
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem 32G
#SBATCH --ntasks=4

module load java/1.8
beagle=/home/sodell/bin/beagle.27Jan18.7e1.jar

#c=$SLURM_ARRAY_TASK_ID
c=10
ref=/group/jrigrp/Share/genotypes/MaizeHapMapV3.2.1/merged_flt_c${c}.vcf.gz
map=/group/jrigrp/Share/annotations/genetic_map/ogut2015/ogut_fifthcM_map_agpv4.txt

java -Xmx32g -jar $beagle nthreads=4 gt=biogemma/Biogemma_WGS_Founders_filtered_snps_diploid_no_tester.vcf.gz out=Biogemma_WGS_Founders_Beagle_Imputed.vcf.gz ref=$ref map=$map chrom=$c impute=true 









