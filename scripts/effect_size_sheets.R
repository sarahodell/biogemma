library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

ses=readRDS('GridLMM/effect_sizes/All_QTL_ES.rds')
snps=ses[which(unlist(lapply(ses,function(x) x$method=="600K_SNP")))]
fps=ses[which(unlist(lapply(ses,function(x) x$method=="Founder_prob")))]
hps=ses[which(unlist(lapply(ses,function(x) x$method=="Haplotype_probs")))]


snp_df=c()
for(i in 1:length(snps)){
  tmp=snps[[i]]
  df=tmp$SE
  df=as.data.frame(df,stringsAsFactors=F)
  df$qtl_id=tmp$qtl_id
  df$snp=tmp$snp
  df$keptsnp=tmp$keptsnp
  df$pos=tmp$pos
  df$method=tmp$method
  df$focal=tmp$focal
  snp_df=rbind(snp_df,df)
}
snp_df=as.data.frame(snp_df,stringsAsFactors=F)
names(snp_df)=c('EffectSize','SE','T-value','Allele','QTL_ID','PeakSNP','KeptSNP','Position','Method','FocalMethod')
fwrite(snp_df,'GridLMM/effect_sizes/SNP_EffectSizes.csv',row.names=F,quote=F,sep=',')

fps_df=c()
for(i in 1:length(fps)){
  tmp=fps[[i]]
  df=tmp$SE[,c('value','se','tvalue','founder')]
  df=as.data.frame(df,stringsAsFactors=F)
  df$qtl_id=tmp$qtl_id
  df$snp=tmp$snp
  df$keptsnp=tmp$keptsnp
  df$pos=tmp$pos
  df$method=tmp$method
  df$focal=tmp$focal
  fps_df=rbind(fps_df,df)
}
fps_df=as.data.frame(fps_df,stringsAsFactors=F)
names(fps_df)=c('EffectSize','SE','T-value','Founder','QTL_ID','PeakSNP','KeptSNP','Position','Method','FocalMethod')
fwrite(fps_df,'GridLMM/effect_sizes/Founder_EffectSizes.csv',row.names=F,quote=F,sep=',')

hps_df=c()
for(i in 1:length(fps)){
  tmp=hps[[i]]
  df=tmp$SE[,c('value','se','tvalue','hapgrp')]
  df=as.data.frame(df,stringsAsFactors=F)
  df$qtl_id=tmp$qtl_id
  df$snp=tmp$snp
  df$keptsnp=tmp$keptsnp
  df$pos=tmp$pos
  df$method=tmp$method
  df$hap_no=tmp$hapgrp
  df$focal=tmp$focal
  hps_df=rbind(hps_df,df)
}
hps_df=as.data.frame(hps_df,stringsAsFactors=F)
names(hps_df)=c('EffectSize','SE','T-value','Haplotype','QTL_ID','PeakSNP','KeptSNP','Position','Method','HaplotypeNumber','FocalMethod')
fwrite(hps_df,'GridLMM/effect_sizes/Haplotype_EffectSizes.csv',row.names=F,quote=F,sep=',')
