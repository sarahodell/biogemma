#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J genoprob
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 800:00:00
#SBATCH --mem 32G
#SBATCH --ntasks=4
#SBATCH --array 0-41%4

module load R/3.3.1

#1 and 9 already completed
sue=( 17 25 33 41 49 57 65 73 81 89 97 105 113 121 129 137 145 153 161 169 177 185 193 201 209 217 225 233 241 249 257 265 273 281 289 297 305 313 321 329 337 )

start=${sue[$SLURM_ARRAY_TASK_ID]}
end=$(($start + 7))
chr=10

echo $start
echo $end
Rscript wgs_genoprobs.R $chr $start $end










