#!/usr/bin/env Rscript

start_marker="AX-91102965"
marker="AX-91102970"
end_marker="AX-91103034"

hap=readRDS('bg8_haplogroup12_reformat.rds')
m=lapply(hap,function(x) x[,c(start_marker,marker)])


saveRDS(m,'AX_91102970_haploprobs.rds')

