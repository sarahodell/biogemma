#!/usr/bin/env Rscript
###Run R/qtl2 on 400 Simulated 16-way MAGIC DH lines 

library('qtl2')

#args=commandArgs(trailingOnly=TRUE)
#c=as.character(args[1])

bg<-read_cross2('MAGICSim_121718_c10.json')

pr <- calc_genoprob(bg,error_prob=0.002,cores=4)
print(dim(pr[[1]]))
saveRDS(pr,"MAGICSim_121718_chr10_genoprobs.rds")