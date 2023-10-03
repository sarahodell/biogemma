#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J magicsim
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=74-97
#SBATCH --mem 8G
#SBATCH --ntasks 1

#Creates a merged vcf files of 400 Simulated MAGIC lines

# Need to rerun 66
module load conda
conda activate /home/sodell/python3.9_env
module load bcftools
module load tabix
module load vcftools

#module load R
#module load python/2.7.14
#
#rep=1
rep=$SLURM_ARRAY_TASK_ID
#founder_path=../../genotypes/600K/individual_founders
#Rscript SimulateMAGIC.R $rep

if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell;
fi

if [ ! -d /scratch/sodell/rep${rep}_tmp ];then
  mkdir /scratch/sodell/rep${rep}_tmp;
fi

#module load conda3
#source activate pandas-env
cd /scratch/sodell/rep${rep}_tmp
founder_path=/home/sodell/projects/biogemma/genotypes/600K/individual_founders
headerfile=/home/sodell/projects/biogemma/genotypes/600K/Biogemma_Founders_600K_Genotypes_AGPv4_no_tester.vcf.gz
build_simvcf=/home/sodell/projects/biogemma/run_magicsim/build_simvcf.py
break_path=/home/sodell/projects/biogemma/run_magicsim/breaktables/MAGIC_DH_Sim_rep${rep}_breaktable_v3.txt
out_path=MAGIC_DHSim_rep${rep}.vcf

echo "Building single line vcfs"

python $build_simvcf $break_path $out_path $headerfile $founder_path --all True
cd /home/sodell/projects/biogemma/run_magicsim
#i=1
#conda deactivate
echo "Ordering and indexing vcfs"
for i in {1..325}; do
    echo Sim${i} > /scratch/sodell/rep${rep}_tmp/Sim${i}.txt
    cat /scratch/sodell/rep${rep}_tmp/Sim${i}_MAGIC_DHSim_rep${rep}.vcf | vcf-sort | bcftools reheader -s /scratch/sodell/rep${rep}_tmp/Sim${i}.txt -o /scratch/sodell/rep${rep}_tmp/Sim${i}_MAGIC_DHSim_rep${rep}_edit.vcf
    bgzip /scratch/sodell/rep${rep}_tmp/Sim${i}_MAGIC_DHSim_rep${rep}_edit.vcf
    tabix -p vcf /scratch/sodell/rep${rep}_tmp/Sim${i}_MAGIC_DHSim_rep${rep}_edit.vcf.gz
    echo /scratch/sodell/rep${rep}_tmp/Sim${i}_MAGIC_DHSim_rep${rep}_edit.vcf.gz >> /scratch/sodell/rep${rep}_tmp/line_names.txt
done
echo "Merging vcfs"
bcftools merge -l /scratch/sodell/rep${rep}_tmp/line_names.txt -m all -Oz > merged_vcfs/MAGIC_DHSimAll_rep${rep}_v3.vcf.gz

#rm *regions.txt
rm -r /scratch/sodell/rep${rep}_tmp
#rm line_names.txt
