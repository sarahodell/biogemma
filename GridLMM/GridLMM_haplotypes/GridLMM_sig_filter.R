#!/usr/bin/env Rscript

## Filtering the p-value results by significance values
args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])

library("data.table")
library("dplyr")


results=fread(sprintf('sig_tables/chr%s_%s_x_%s.txt',chr,pheno,env))

k = data.frame(results %>% group_by(hapgrp) %>% summarize(size=length(bin)))
cutoff=-log10(0.05/sum(k$size))

sig=c()
for(x in 1:15){
  sig=rbind(sig,results[(results$hapgrp==k$hapgrp[x]) & (results$log10p>=cutoff),])
}

if(dim(sig)[1]!=0){
  print(dim(sig))
  fwrite(sig,sprintf('sig_tables/chr%s_%s_x_%s_bonf_cutoff.txt',chr,pheno,env),quote=F,row.names=F,sep=',')
}