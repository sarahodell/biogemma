#!/bin/bash

chr=$1

module load bedtools

if [ -f ibd_segments/comparison/c${chr}_germline_overlap.bed ]; then rm ibd_segments/comparison/c${chr}_germline_overlap.bed; fi
if [ -f ibd_segments/comparison/c${chr}_600K_ibdseq_overlap.bed ]; then rm ibd_segments/comparison/c${chr}_600K_ibdseq_overlap.bed; fi
if [ -f ibd_segments/comparison/c${chr}_germline_wgs_ibdseq_overlap.bed ]; then rm ibd_segments/comparison/c${chr}_germline_wgs_ibdseq_overlap.bed; fi
if [ -f ibd_segments/comparison/c${chr}_germline_600K_RefinedIBD_overlap.bed ]; then rm ibd_segments/comparison/c${chr}_gerline_600K_RefinedIBD_overlap.bed; fi

echo "GERMLINE 600K vs. GERMLINE WGS"
while read pair; do
  bedtools intersect -a ibd_segments/comparison/tmp/${pair}_wgs_c${chr}.bed -b ibd_segments/comparison/tmp/${pair}_600K_c${chr}.bed >> ibd_segments/comparison/c${chr}_germline_overlap.bed
done < ibd_segments/comparison/c${chr}_germline_pair_list.txt

echo "GERMLINE 600K vs. IBDSeq WGS"
while read pair; do
  bedtools intersect -a ibd_segments/comparison/tmp/${pair}_600K_c${chr}.bed -b ibd_segments/comparison/tmp/${pair}_wgs_ibdseq_c${chr}.bed >> ibd_segments/comparison/c${chr}_600K_ibdseq_overlap.bed
done < ibd_segments/comparison/c${chr}_600K_ibdseq_pair_list.txt

echo "GERMLINE WGS vs. IBDSeq WGS"
while read pair; do
  bedtools intersect -a ibd_segments/comparison/tmp/${pair}_wgs_c${chr}.bed -b ibd_segments/comparison/tmp/${pair}_wgs_ibdseq_c${chr}.bed >> ibd_segments/comparison/c${chr}_germline_wgs_ibdseq_overlap.bed
done < ibd_segments/comparison/c${chr}_germline_wgs_ibdseq_pair_list.txt

echo "GERMLINE 600K vs. RefinedIBD"

while read pair; do
  bedtools intersect -a ibd_segments/comparison/tmp/${pair}_600K_c${chr}.bed -b ibd_segments/comparison/tmp/${pair}_RefinedIBD_c${chr}.bed >> ibd_segments/comparison/c${chr}_germline_600K_RefinedIBD_overlap.bed
done < ibd_segments/comparison/c${chr}_germline_600K_RefinedIBD_pair_list.txt
