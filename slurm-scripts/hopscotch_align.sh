#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/genotypes/vgt1
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%j-out.txt
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%j-error.txt
#SBATCH -J mite
#SBATCH -t 96:00:00
#SBATCH --mem 64G
#SBATCH --ntasks=8
#SBATCH --array=1-17%3

declare -a sample_list=("blank" "A632_usa" "A654_inra" "B73_inra" "B96" "C103" "CO255_inra" "D105_inra" "DK63" "EP1_inra" "F492" "FV252_inra" "FV2_inra" "MBS847" "ND245" "OH43_inra" "VA85" "W117_inra" )

SAMPLE=$(echo ${sample_list[$SLURM_ARRAY_TASK_ID]})
#SAMPLE=A632_usa
echo $SAMPLE
module load samtools
module load bamtools

FASTA=sequences/C22_4_mite.fa
EMPTYFA=sequences/N28_no_mite.fa
CPU=8

#SAMPLE=$SAMPLES[$SLURM_ARRAY_TASK_ID]

READ1=/scratch/sodell/Sample_${SAMPLE}_R1_fastq.txt.gz
READ2=/scratch/sodell/Sample_${SAMPLE}_R2_fastq.txt.gz

ls /scratch/sodell/*

#samtools faidx $FASTA
#samtools faidx $EMPTYFA
bwa index $FASTA
bwa index $EMPTYFA

bwa mem -t $CPU -a $EMPTYFA $READ1 $READ2 | samtools view -Su - | bamtools filter -script filter_mapped_and_pairs.bamtools.json | samtools sort - > bams/vgt1/${SAMPLE}.no_hopscotch.sorted.bam
samtools index bams/vgt1/${SAMPLE}.no_hopscotch.sorted.bam


#if [ -f /home/sodell/projects/biogemma/genotypes/vgt1/bams/vgt1/${SAMPLE}.no_hopscotch.sorted.bam ]
#then
samtools sort -n bams/vgt1/${SAMPLE}.no_hopscotch.sorted.bam | samtools fastq - | bwa mem -t $CPU -p $FASTA - | samtools sort - > bams/vgt1/${SAMPLE}.with_hopscotch.sorted.bam

#fi
samtools index bams/vgt1/${SAMPLE}.with_hopscotch.sorted.bam
#samtools index bams/vgt1/${SAMPLE}.no_hopscotch.sorted.bam

#Grab reads with proper pairs
#samtools view -b -F 12 bams/vgt1/${SAMPLE}.no_hopscotch.sorted.bam -o bams/vgt1/${SAMPLE}_no_hopscotch_paired_only.bam
#samtools index bams/vgt1/${SAMPLE}_no_hopscotch_paired_only.bam


#Remove reads that don't span the MITE

samtools view -b -F 12 bams/vgt1/${SAMPLE}.with_hopscotch.sorted.bam -o bams/vgt1/${SAMPLE}_with_hopscotch_paired_only.bam
samtools index bams/vgt1/${SAMPLE}_with_hopscotch_paired_only.bam

#Count reads that map to the MITE, compare based off alignment to emptyfa and fasta
