#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J hmptovcf
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load GATK/3.6

java -jar /share/apps/GATK-3.6/GenomeAnalysisTK.jar -T VariantsToVCF -R /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa -o biogemma/Biogemma_Founders_600K_AGPv4_ref.vcf --variant:RawHapMap biogemma/Biogemma_Founders.hmp.txt 
