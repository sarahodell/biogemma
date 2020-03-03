#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

library('data.table')
library('dplyr')
library('ggplot2')

all_reps=c()
for(i in 1:10){
  r=readRDS(sprintf('test_models/chr%.0f_%s_x_%s_600KSNP_1000rep.rds',i,pheno,env),data.table=F)
  df=c()
  df=sapply(seq(1,1000),function(x) rbind(df,unlist(r[[x]])))
  df=t(df)
  df=as.data.frame(df)
  names(df)=c('chr','replicate','pval')
  df=df[!is.na(df$pval),]
  all_reps=rbind(all_reps,df)
}

minp = all_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
minp=as.data.frame(minp)

threshold=quantile(minp$pval,0.05,lower.tail=T)
print(threshold)
print(-log10(threshold))

method="600K_SNP"

line=data.table(phenotype=pheno,environment=env,method=method,threshold=-log10(threshold))
fwrite(line,file="../../threshold_table.txt",sep=',',append=T)


