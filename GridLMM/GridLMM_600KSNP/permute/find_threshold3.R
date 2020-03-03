#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

library('data.table')
library('dplyr')
library('ggplot2')

options(scipen=20)

all_reps=c()
for(c in 1:10){
  chr_reps=c()
  reps=seq(1,1000)
  for(j in 1:10){
    r=readRDS(sprintf('test_models/chr%.0f_%s_x_%s_600KSNP_%.0frep.rds',c,pheno,env,j))
    if(sum(sapply(seq(1,100),function(x) is.null(r[[x]]))) != 0){print("Null Error")}
    df=c()
    df=sapply(seq(1,100),function(x) rbind(df,unlist(r[[x]])))
    df=t(df)
    df=as.data.frame(df)
    names(df)=c('chr','replicate','pval')
    df=df[!is.na(df$pval),]
    chr_reps=rbind(chr_reps,df)
  }
  chr_reps$chr=c
  chr_reps$replicate=reps
  all_reps=rbind(all_reps,chr_reps)
}

minp = all_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
minp=as.data.frame(minp)

threshold=quantile(minp$pval,0.05,lower.tail=T)
print(threshold)
print(-log10(threshold))

method="600K_SNP"

line=data.table(phenotype=pheno,environment=env,method=method,threshold=-log10(threshold))
fwrite(line,file="../../threshold_table.txt",sep=',',append=T)


#png(sprintf('%s_x_%s_perm_1000_pval_dist.png',pheno,env))
#print(ggplot(minp,aes(x=pval)) + geom_histogram() + geom_vline(xintercept=threshold))
#dev.off()

#png(sprintf('%s_x_%s_perm_1000_log10pval_dist.png',pheno,env))
#print(ggplot(minp,aes(x=-log10(pval))) + geom_histogram() + geom_vline(xintercept=-log10(threshold)))
#dev.off()





