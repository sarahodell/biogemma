#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J mpileup
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=7

export PATH=/home/sodell/bin/bin/:$PATH 

#while read line;do
#    echo -e /group/jrigrp/Share/sequences/biogemma_magic/Sample_${line}.sorted.bam >> bamfiles.txt
#done<biogemma/founders/biogemma_founders.txt

#samtools mpileup -uBI -l biogemma/Axiom600K_v4.bed -f /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa -b bamfiles.txt | bcftools call -c - > BiogemmaFoundersWGS_to_Axiom.vcf 
     
#bcftools mpileup -Ou -I -Q 0 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP -b bamfiles.txt -f /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa -R biogemma/axiom600K_positions.txt | bcftools call -Oz -m > BiogemmaFoundersWGS_to_Axiom_2.vcf.gz

bcftools mpileup -Ou -I -Q 0 -a INFO/AD,INFO/ADF,INFO/ADR,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP -f /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa -b bamfiles.txt -T biogemma/axiom600K_positions.txt | bcftools call -Oz -m > BiogemmaFoundersWGS_to_Axiom_2.vcf.gz


