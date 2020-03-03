#!/bin/bash

for c in {1..10};do
    by=43
    start=2
    file=Biogemma_080618/Biogemma_WGS_geno_chr${c}.csv
    echo $file
    for b in {1..8};do
	end=$(($by * $b))
	#echo $start
	#echo $end
	sed -n -e 1p -e ${start},${end}p $file > Biogemma_080618/Biogemma_WGS_geno_chr${c}_${b}.csv  
	start=$(($end + 1))
    done
done
