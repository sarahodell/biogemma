#!/bin/bash -l
#SBATCH -D /home/sodell/projects/impute
#SBATCH -J assessproja
#SBATCH -o /home/sodell/projects/impute/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/impute/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

outfile='Predicted_BiogemmaMagic_4.txt'
proja='MAGICSim_062018/FILLIN/MAGICSim_Biogemma_062018_Imputed4.pa.txt'
actual='MAGICSim_062018/Biogemma_MAGICSim.txt'

python pscripts/assess_proja.py $proja $outfile
python pscripts/intersect.py $actual $outfile












