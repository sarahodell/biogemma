#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
#chr=as.character(args[[3]])
#cores=as.numeric(args[[4]])

#date=format(Sys.time(),'%m%d%y')

library('ggplot2')
library('data.table')
library('dplyr')

base_list=c(10,9,9,10,7,10,6,9,9,10)
#base_list=c(7,10,6,7,9,9,8,9,8,2)
all_chroms=c()

for(i in 1:10){
   pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
   ID=c()
   pos=c()
   bin=c()
   pvalues=c()
   hapgrp=c()
   base=base_list[i]
   for(h in base:16){
      mod=readRDS(sprintf('GridLMM_haplotypes/models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',i,h,pheno,env))
      p=match(mod$X_ID,pmap$marker)
      phy_pos=pmap[p,]$pos
      ID=c(ID,mod$X_ID)
      pos=c(pos,phy_pos)
      pvalues=c(pvalues,mod$p_value_ML)
      hapgrp=c(hapgrp,rep(h,length(mod$p_value_ML)))
      bin=c(bin,seq(1,length(mod$p_value_ML)))
   }
   gwas=data.frame(chr=i,hapgrp=hapgrp,bin=bin,ID=ID,pos=pos,pvalues=pvalues,stringsAsFactors=F)
   all_chroms=rbind(all_chroms,gwas)
}
all_chroms$chr=as.factor(all_chroms$chr)
all_chroms$bin=as.factor(all_chroms$bin)
all_chroms$hapgrp=as.factor(all_chroms$hapgrp)
all_chroms$log10p = -log10(all_chroms$pvalues)

size=dim(all_chroms)[1]

thresh_table=fread('threshold_table.txt',data.table=F,stringsAsFactors=F)
rec=thresh_table$phenotype==pheno & thresh_table$environment==env & thresh_table$method=="haplotype_probs"
cutoff=thresh_table[rec,]$threshold
print(cutoff)

all_chroms$sig = all_chroms$log10p>= cutoff
label=sprintf('5%% Permutation Significance Threshold = %.3f',cutoff)

png(sprintf('GridLMM_haplotypes/images/%s_x_%s_manhattan_sig.png',pheno,env),width=960,height=680)
theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(size=26),axis.title=element_text(size=14,face="bold"))
theme_update(panel.background=element_blank())

print(ggplot(all_chroms,aes(x=pos/1e6,y=log10p)) + geom_point(aes(color=sig)) + scale_color_manual(breaks=all_chroms$sig,values=c("FALSE"="black","TRUE"="red")) + facet_grid(.~chr,scales="free") + ggtitle(sprintf("%s in %s Using Haplotype Probabilities",pheno,env)) + xlab("Chromosome (Mb)") + ylab("-log10(P-Value)") + guides(color=F) + labs(caption=label))
dev.off()
