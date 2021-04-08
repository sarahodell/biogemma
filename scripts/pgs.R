#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
env=as.character(args[[1]])
cores=as.numeric(args[[2]])

library('data.table')
#library('tidyverse')
library('ggplot2')
library('parallel')
library('MASS')

c=10

base_list=c(7,7,8,6,7,6,7,7,6,8)

K=fread('GridLMM/K_matrices/K_matrix_chr10.txt',data.table=F)
inds=K$V1
s_pgs=rep(0,length(inds))
f_pgs=rep(0,length(inds))
h_pgs=rep(0,length(inds))

geno_list=vector("list",length=10)
ses_list=vector("list",length=10)
for(c in 1:10){
  #SNPs
  geno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%.0f_binary.csv',c),data.table=F)
  ses=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%.0f_male_flowering_d6_x_%s_600KSNP_ML.rds',c,env))
  ses=ses$results
  betas = ses[,grep('beta',colnames(ses))]
  betas$beta.2 = betas$beta.2+betas$beta.1
  betas[betas$beta.1 < 0.05 & betas$beta.1 > -0.05,]$beta.1 = 0
  betas[betas$beta.2 < 0.05 & betas$beta.2 > -0.05,]$beta.2 = 0
  rownames(betas)=ses$X_ID
  rownames(geno)=geno$ind
  geno=geno[inds,ses$X_ID]
  geno_list[[c]]=geno
  ses_list[[c]]=betas
}
  #Founders
  #fes=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_male_flowering_d6_x_%s_founderprobs.rds',c,env))

  #Haplotype
  #base=base_list[c]
  #for(h in base:16){
  #  hapfile=sprintf('GridLMM/GridLMM_haplotypeprobs/models/Biogemma_chr%.0f_haplogrp%.0f_male_flowering_d6_x_%s.rds',c,h,env)
  #  if file.exists(hapfile){
  #    hes=readRDS(hapfile)
  #  }
  #}


chr_pgs<-function(c,ind){
  geno=geno_list[[c]]
  bts=ses_list[[c]]
  nmarkers=dim(geno)[2]
  t=geno[ind,]
  b=paste0('beta.',t+1)
  grab_row=rownames(bts)
  grab_col=b
  pg=mapply(function(i,j) bts[i,j], grab_row, grab_col)
  #pg=sapply(seq(1,length(t)), function(x) betas[names(t)[x],b[x]])
  #pg=sapply(seq(1,nmarkers), function(i) betas[rownames(betas)[i],paste0('beta.',(geno[ind,)+1)])
  pgsum=sum(pg)
  names(pgsum)=ind
  return(pgsum)
}
print(system.time({
  results=chr_pgs(10,'EB.09S.H.00002')
}))


calc_pgs<-function(ind){
  chr_res=0
  chroms=seq(1,10)
  for(c in 1:10){
    chr_res=chr_res+chr_pgs(c,ind)
  }
  #ses=ses[match(ses$X_ID),]
  #s_pgs=s_pgs + snp_pg
  names(chr_res)=ind
  return(chr_res)
}
#print(system.time({
#  results=calc_pgs('EB.09S.H.00002')
#}))

n_inds=inds

print(system.time({
results=mclapply(n_inds,calc_pgs,mc.cores=cores)
}))

all_pgs=data.frame(inds=inds,s_pgs=unlist(results),stringsAsFactors=F)
fwrite(all_pgs,'GridLMM/pgs/SNP_dta_polygenic_scores.txt',quote=F,sep='\t',row.names=F)
