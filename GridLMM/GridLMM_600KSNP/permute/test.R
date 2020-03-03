#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])

r=readRDS(sprintf('test_models/chr%s_%s_x_%s_600KSNP_1000rep.rds',chr,pheno,env))

result=sum(sapply(seq(1,1000),function(x) is.null(r[[x]])))

if(result!=0){
  print(result)
}