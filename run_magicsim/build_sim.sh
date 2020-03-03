#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute/for_asher
#SBATCH -J magicsim
#SBATCH -o /home/sodell/projects/impute/for_asher/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/for_asher/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --mem 32G
#SBATCH --ntasks 4

#Creates a merged vcf files of 400 Simulated MAGIC lines

module load bcftools
module load tabix
module load vcftools
module load R/3.5.2
module load python/2.7.14 
module load conda3/1.0

source activate pandas-env
#Whatever today's date is
today="012519"

#Rscript SimulateMAGIC.R

mkdir tmp

python build_simvcf.py MAGIC_DHSim_${today}_breakpoints.txt MAGIC_DHSim_${today}.vcf founders --all True

for i in {1..400}; do
    echo Sim${i} > tmp/Sim${i}.txt
    cat tmp/Sim${i}_MAGIC_DHSim_${today}.vcf | vcf-sort | bcftools reheader -s tmp/Sim${i}.txt -o tmp/Sim${i}_MAGIC_DHSim_${today}_edit.vcf
    bgzip tmp/Sim${i}_MAGIC_DHSim_${today}_edit.vcf
    tabix -p vcf tmp/Sim${i}_MAGIC_DHSim_${today}_edit.vcf.gz
    echo tmp/Sim${i}_MAGIC_DHSim_${today}_edit.vcf.gz >> line_names.txt
done

bcftools merge -l line_names.txt -m all -Oz > MAGIC_DHSimAll_${today}_chr10.vcf.gz

rm *regions.txt
rm -r tmp
rm line_names.txt
