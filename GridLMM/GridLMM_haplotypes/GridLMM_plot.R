#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])

#date=format(Sys.time(),'%m%d%y')

library('ggplot2')
library('data.table')
library('dplyr')

bin=c()
pvalues=c()
hapgrp=c()
for(i in 2:16){
    mod=readRDS(sprintf('models/Biogemma_chr%s_haplogrp%.0f_%s_x_%s.rds',chr,i,pheno,env))
    pvalues=c(pvalues,mod$p_value_ML)
    hapgrp=c(hapgrp,rep(i,length(mod$p_value_ML)))
    bin=c(bin,seq(1,length(mod$p_value_ML)))
}
gwas=data.frame(hapgrp=hapgrp,bin=bin,pvalues=pvalues)
gwas$bin=as.factor(gwas$bin)
gwas$hapgrp=as.factor(gwas$hapgrp)
gwas$log10p = -log10(gwas$pvalues)

#k=as.data.frame(gwas %>% group_by(hapgrp) %>% summarize(size=length(bin)))
size=dim(gwas)[1]

fwrite(gwas,sprintf('sig_tables/chr%s_%s_x_%s.txt',chr,pheno,env),quote=F,row.names=F,sep='\t')

png(sprintf('images/chr%s_%s_x_%s_manhattan.png',chr,pheno,env),width=960,height=680)
print(ggplot(gwas,aes(x=bin,y=log10p)) + geom_point(aes(color=hapgrp)) + geom_hline(aes(yintercept=-log10(0.05/size)),color="red") + facet_grid(.~hapgrp,scales="free") + ggtitle(sprintf("%s in %s on Chromosome %s",pheno,env,chr)) + xlab("Haplotype Group") + ylab("-log10(p-value)") + theme_classic())
dev.off()



