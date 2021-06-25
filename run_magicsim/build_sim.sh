#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J magicsim
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-100
#SBATCH --mem 3G
#SBATCH --ntasks 1

#Creates a merged vcf files of 400 Simulated MAGIC lines

module load bcftools
module load tabix
module load vcftools
module load R
module load python/2.7.14
#module load conda3

#source activate pandas-env
#Whatever today's date is
#today="012519"
#rep=
rep=$SLURM_ARRAY_TASK_ID

Rscript SimulateMAGIC.R $rep

#mkdir tmp

#python build_simvcf.py MAGIC_DHSim_breakpoints.txt MAGIC_DHSim.vcf founders --all True

#for i in {1..400}; do
#    echo Sim${i} > tmp/Sim${i}.txt
#    cat tmp/Sim${i}_MAGIC_DHSim.vcf | vcf-sort | bcftools reheader -s tmp/Sim${i}.txt -o tmp/Sim${i}_MAGIC_DHSim_${today}_edit.vcf
#    bgzip tmp/Sim${i}_MAGIC_DHSim_edit.vcf
#    tabix -p vcf tmp/Sim${i}_MAGIC_DHSim_edit.vcf.gz
#    echo tmp/Sim${i}_MAGIC_DHSim_edit.vcf.gz >> line_names.txt
#done

#bcftools merge -l line_names.txt -m all -Oz > MAGIC_DHSimAll_chr10.vcf.gz

#rm *regions.txt
#rm -r tmp
#rm line_names.txt
