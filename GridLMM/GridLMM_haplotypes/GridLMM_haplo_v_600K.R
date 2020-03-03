#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

library('ggplot2')
library('data.table')
library('dplyr')

options(scipen=99)
comp=c()

rsq <- function(x,y){
   return(cor(x,y) ^2)
}

thresh_table=fread('threshold_table.txt',data.table=F,stringsAsFactors=F)
rec1=thresh_table$phenotype==pheno & thresh_table$environment==env & thresh_table$method=="haplotype_probs"
  rec2=thresh_table$phenotype==pheno & thresh_table$environment==env & thresh_table$method=="600K_SNP"
hap_cutoff=thresh_table[rec1,]$threshold
snp_cutoff=thresh_table[rec2,]$threshold

pdf(sprintf('method_comparison/%s_x_%s_haplo_x_600K_comparison.pdf',pheno,env))
theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(size=14),axis.title=element_text(size=12,face="bold"))
theme_update(panel.background=element_blank())
for(i in 1:10){
   pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
   ID1=c()
   pos1=c()
   pvalues1=c()
   hapgrp=c()
   for(h in 2:16){
      mod=readRDS(sprintf('models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',i,h,pheno,env))
      p=match(mod$X_ID,pmap$marker)
      phy_pos=pmap[p,]$pos
      ID1=c(ID1,mod$X_ID)
      pos1=c(pos1,phy_pos)
      pvalues1=c(pvalues1,mod$p_value_ML)
      hapgrp=c(hapgrp,rep(h,length(mod$p_value_ML)))
   }
   hap_gwas=data.frame(chr=i,hapgrp=hapgrp,ID=ID1,pos=pos1,pvalues=pvalues1,stringsAsFactors=F)
   hap_gwas=hap_gwas[!is.na(hap_gwas$pvalues),]
   hap_gwas=hap_gwas[!is.na(hap_gwas$pos),]
   hap_gwas$log10p=-log10(hap_gwas$pvalues)
   mod2=readRDS(sprintf('GridLMM_600KSNP/models/chr%.0f_%s_x_%s_600KSNP.rds',i,pheno,env))
   p2=match(mod2$results$X_ID,pmap$marker)
   phy_pos=pmap[p2,]$pos
   ID2=mod2$results$X_ID
   snp_gwas=data.frame(chr=i,ID=ID2,pos=phy_pos,pvalues=mod2$results$p_value_REML,stringsAsFactors=F)
   snp_gwas=snp_gwas[!is.na(snp_gwas$pvalues),]
   snp_gwas=snp_gwas[!is.na(snp_gwas$pos),]
   snp_gwas$log10p=-log10(snp_gwas$pvalues)
   
   max_pos=max(pmap$pos)
   chr_comp=c()
   start=0
   end=2e6
   bin=1
   for(k in seq(1,as.integer(max_pos/2e6))){
      snp_bin=max(snp_gwas[snp_gwas$pos>start & snp_gwas$pos <=end,]$log10p)
      hap_bin=max(hap_gwas[hap_gwas$pos>start & hap_gwas$pos <=end,]$log10p)
      if((is.na(snp_bin)==F) & (is.infinite(snp_bin)==F) & (is.na(hap_bin)==F) & (is.infinite(hap_bin)==F)){
         chr_comp=rbind(chr_comp,data.frame(bin=sprintf('%s_%s',i,bin),bin_start=start,bin_end=end,chr=i,snp_minp=snp_bin,hap_minp=hap_bin,stringsAsFactors=F))
      }
      start=end+1
      end=end+2e6
      bin=bin+1
   }

   r2=rsq(chr_comp$snp_minp,chr_comp$hap_minp)
   label=sprintf("R-Squared = %.3f, Bin size of 2Mb",r2)

   print(ggplot(chr_comp,aes( x=snp_minp,y=hap_minp)) + geom_point() + geom_vline(xintercept=snp_cutoff) + geom_hline(yintercept=hap_cutoff) + ggtitle(sprintf("%s in %s Haplotype vs 600K p-values, Chr %.0f",pheno,env,i)) + xlab("600K Max -log10(p-value) per Bin") + ylab("Haplotype Max -log10(p-value) per Bin") + guides(color=F) + labs(caption=label))
  comp=rbind(comp,chr_comp)
}
dev.off()

fwrite(comp,sprintf('method_comparison/%s_x_%s_haplo_x_600KSNP_2MBbins.txt',pheno,env),quote=F,sep='\t',row.names=F)


