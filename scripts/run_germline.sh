#!/bin/bash

germline=/home/sodell//bin/germline-1-5-3/./germline 

for i in {1..10}; do
    $germline -input Biogemma_Founders_600K_chr${i}.ped Biogemma_Founders_600K_chr${i}.map -min_m 1 -err_hom 4 -w_extend -output germline/Biogemma_Founders_germline_IBD_chr${i} 
done

#i=10
#$germline -input Biogemma_Founders_600K_chr${i}.ped Biogemma_Founders_600K_chr${i}.map -min_m 1 -err_hom 4 -w_extend -output germline/Biogemma_Founders_germline_IBD_chr${i}

