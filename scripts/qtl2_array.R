#!/usr/bin/env Rscript
###Run R/qtl2 on 344 16-way MAGIC DH lines 

library('qtl2')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])
cores=as.character(args[[2]])

control_file=sprintf('Biogemma_c%s.json',c)
outfile=sprintf("bg%s_genoprobs_010319.rds",c)

bg<-read_cross2(control_file)

pr <- calc_genoprob(bg,error_prob=0.002,cores=cores)
#print(dim(pr[[1]]))

saveRDS(pr,outfile)