#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/run_magicsim
#SBATCH -J plink
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 12:00:00
#SBATCH --array=1-100
#SBATCH --ntasks=16
#SBATCH --mem 127G

module load plink

if [ ! -d /scratch/sodell ]; then
  mkdir /scratch/sodell;
fi
rep=$SLURM_ARRAY_TASK_ID
plink --vcf merged_vcfs/MAGIC_DHSimAll_rep${rep}.vcf.gz --out plinkfiles/MAGIC_DHSimAll_rep${rep}
#rep=1

plink --threads 16 --bfile plinkfiles/MAGIC_DHSimAll_rep${rep} --r2 with-freqs inter-chr --ld-window-r2 0.9 --out /scratch/sodell/MAGIC_DHSimAll_rep${rep}_r2_interchrom

if [ ! -f inter_ld_count.txt ]; then
  touch inter_ld_count.txt;
fi

count=$(awk '$1!=$5{print $1,$5}' /scratch/sodell/MAGIC_DHSimAll_rep${rep}_r2_interchrom.ld | wc -l)
count2="$(($count-1))"

echo $rep $count2 >> inter_ld_count.txt

rm /scratch/sodell/MAGIC_DHSimAll_rep${rep}_r2_interchrom.ld
