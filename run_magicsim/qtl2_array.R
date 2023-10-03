#!/usr/bin/env Rscript
###Run R/qtl2 on 344 16-way MAGIC DH lines

library('broman')
library('qtl2')
library('qtlcharts')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])
rep=as.character(args[2])
cores=as.numeric(args[[3]])

control_file=sprintf('MAGIC_DHSim_rep%s_c%s_v3.json',rep,c)
outfile=sprintf("MAGIC_DHSim_rep%s_c%s_genoprobs_v3.rds",rep,c)

bg<-read_cross2(control_file)
bg<-drop_nullmarkers(bg)

prcl <- calc_genoprob(bg,error_prob=0.002,cores=cores)
pr_clean <- clean_genoprob(prcl)
#print(dim(pr[[1]]))

saveRDS(pr_clean,outfile)
