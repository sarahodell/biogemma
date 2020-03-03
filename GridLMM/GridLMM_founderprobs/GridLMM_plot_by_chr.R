#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

library('ggplot2')
library('cowplot')
library('data.table')
library('dplyr')

all_chroms=c()

for(i in 1:10){
   pmap=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F,stringsAsFactors=F)
   mod=readRDS(sprintf('models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',i,pheno,env))
   p=match(mod$X_ID,pmap$marker)
   phy_pos=pmap[p,]$pos
   ID=mod$X_ID
   gwas=data.frame(chr=i,pos=phy_pos,pvalues=mod$p_value_ML,stringsAsFactors=F)
   all_chroms=rbind(all_chroms,gwas)
}
all_chroms$chr=as.factor(all_chroms$chr)
all_chroms$log10p = -log10(all_chroms$pvalues)
size=dim(all_chroms)[1]
all_chroms=all_chroms[order(all_chroms$pvalues),]
rownames(all_chroms)=seq(1,size)


thresh_table=fread('../threshold_table.txt',data.table=F,stringsAsFactors=F)
rec=thresh_table$phenotype==pheno & thresh_table$environment==env & thresh_table$method=="founder_probs"
cutoff=thresh_table[rec,]$threshold
print(cutoff)

all_chroms$sig=all_chroms$log10p >= cutoff
label<-sprintf("5%% Permutation Significance Threshold = %.3f",cutoff)


png(sprintf('images/%s_x_%s_manhattan.png',pheno,env),width=960,height=680)
theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(size=26),axis.title=element_text(size=14,face="bold"))
theme_update(panel.background=element_blank())
print(ggplot(all_chroms,aes(x=pos/1e6,y=log10p)) + geom_point(aes(color=sig)) + scale_color_manual(breaks=all_chroms$sig,values=c("FALSE"="black","TRUE"="red")) + facet_grid(.~chr,scales="free") + ggtitle(sprintf("%s in %s Using Founder Probabiliites",pheno,env)) + xlab("Chromosome (Mb)") + ylab("-log10(P-Value)") + guides(color=F) + labs(caption=label))
dev.off()



