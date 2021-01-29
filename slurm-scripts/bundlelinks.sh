#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma
#SBATCH -J bundlelinks
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=8
#SBATCH --mem 60G


module load perlbrew
module load circos

circostools=/home/sodell/bin/circos-tools-0.23/tools/bundlelinks
linkfile=stats/ld_decay/circos/ld_links_interchrom.txt
outfile=stats/ld_decay/circos/ld_bundled_links.txt

cd $circostools

#bin/bundlelinks -links $linkfile > $outfile

./run
