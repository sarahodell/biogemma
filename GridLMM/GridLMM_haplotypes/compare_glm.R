#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
c=as.character(args[[2]])

library('data.table')



envs=c("ALL","SZEGED_2017_OPT","STPAUL_2017_WD","NERAC_2016_WD","GRANEROS_2015_OPT","BLOIS_2017_OPT","BLOIS_2014_OPT")
if(pheno=="male_flowering_days" | pheno=="female_flowering_days"){
  envs=c("ALL","NERAC_2016_WD","GRANEROS_2015_OPT","BLOIS_2014_OPT")
}
### Read in GridLMM model output for the three methods and write out to a single file
pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F,stringsAsFactors=F)

thresh_table=fread('threshold_table.txt',data.table=F,stringAsFactors=F)

find_threshold<function(pheno,env,method){
   cond=thresh_table$phenotype==pheno & thresh_table$environment==env & thresh_table$method==method
   threshold=thresh_table[cond,]$cutoff
   return(threshold)
}

print(c)

all_results=c()

# 600K SNP
for(e in envs){
  threshold=find_threshold(pheno,e,'600K_SNP')
  mod=readRDS(sprintf('GridLMM_600KSNP/models/chr%s_%s_x_%s_600KSNP.rds',c,pheno,e))$results
  mod=mod[!is.na(mod$p_value_REML),]
  submod=data.frame(marker=mod$X_ID,p_value=mod$p_value_REML,chr=c,method="600K_SNP",sig=-log10(mod$p_value_REML)>=threshold,env=e,stringsAsFactors=F)
#  print(dim(submod))
  all_results=rbind(all_results,submod)
}
# Founder probabilities
for(e in envs){
  threshold=find_threshold(pheno,e,'founder_probs')
  mod=readRDS(sprintf('GridLMM_founderprobs/models/Biogemma_chr%s_%s_x_%s_founderprobs_adjusted.rds',c,pheno,e))
  mod=mod[!is.na(mod$p_value_ML),]
  submod=data.frame(marker=mod$X_ID,p_value=mod$p_value_ML,chr=c,method="founder_probs",sig=-log10(mod$p_value_ML)>=threshold,env=e,stringsAsFactors=F)
#  print(dim(submod))
  all_results=rbind(all_results,submod)
}


# Haplotype probabilities
for(e in envs){
  threshold=find_threshold(pheno,e,'haplotype_probs')
  for(h in 2:16){
    mod=readRDS(sprintf('models/Biogemma_chr%s_haplogrp%.0f_%s_x_%s_adjusted.rds',c,h,pheno,e))
    mod=mod[!is.na(mod$p_value_ML),]
    submod=data.frame(marker=mod$X_ID,p_value=mod$p_value_ML,chr=c,method="haplotype_probs",sig=-log10(mod$p_value_ML)>=threshold,env=e,stringsAsFactors=F)
    all_results=rbind(all_results,submod)
  }
}

all_results=as.data.frame(all_results)
names(all_results)=c('marker','p_value','chr','method','sig','env')

final_results=merge(all_results,pmap,by.x='marker',by.y='marker',all.x=T)

fwrite(final_results,sprintf('method_comparison/chr%s_%s_all_methods.txt',c,pheno),row.names=F,quote=F)

sig_only=final_results[final_results$sig==T,]
sig_only=sig_only[order(sig_only$pos),]
rownames(sig_only)=seq(1,dim(sig_only)[1])

fwrite(sig_only,sprintf('method_comparison/chr%s_%s_all_methods_sig_only.txt',c,pheno),row.names=F,quote=F)