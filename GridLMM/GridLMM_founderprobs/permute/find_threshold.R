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
   r=fread(sprintf('test_models/chr%.0f_%s_x_%s_founderprobs_1000rep_max_pvalues.txt',c,pheno,env),data.table=F)
   r=r[!is.na(r$pval),]
   tmp=data.frame(chr=c,replicate=r$replicate,pval=r$pval,stringsAsFactors=F)
   all_reps=rbind(all_reps,tmp)
}

fwrite(all_reps,sprintf('max_reps/%s_x_%s_rep1000_max_pvalues.txt',pheno,env),quote=F,row.names=F,sep='\t')

minp = all_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
minp=as.data.frame(minp)

threshold=quantile(minp$pval,0.05,lower.tail=T)
print(threshold)
print(-log10(threshold))

png(sprintf('%s_x_%s_perm_1000_pval_dist.png',pheno,env))
print(ggplot(minp,aes(x=pval)) + geom_histogram() + geom_vline(xintercept=threshold))
dev.off()

png(sprintf('%s_x_%s_perm_1000_log10pval_dist.png',pheno,env))
print(ggplot(minp,aes(x=-log10(pval))) + geom_histogram() + geom_vline(xintercept=-log10(threshold)))
dev.off()





