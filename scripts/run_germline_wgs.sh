#!/bin/bash

germline=/home/sodell/bin/germline-1-5-3/./germline 


for i in {1..10}; do
    input_map=/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr${i}.map
    input_ped=/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr${i}.ped
    $germline -input $input_ped $input_map -w_extend -min_m 1 -err_hom 4 -output ibd_segments/germline/WGS/Biogemma_Founders_WGS_germline_IBD_chr${i} 
done

#i=10
#input_map=/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr${i}.map
#input_ped=/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr${i}.ped
#$germline -input $input_ped $input_map -min_m 1 -err_hom 4 -w_extend -output ibd_segments/germline/WGS/Biogemma_Founders_WGS_germline_IBD_chr${i}

