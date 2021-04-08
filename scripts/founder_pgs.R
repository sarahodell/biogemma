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

#base_list=c(7,7,8,6,7,6,7,7,6,8)

K=fread('GridLMM/K_matrices/K_matrix_chr10.txt',data.table=F)
inds=K$V1
#s_pgs=rep(0,length(inds))
#f_pgs=rep(0,length(inds))
#h_pgs=rep(0,length(inds))

geno_list=vector("list",length=10)
fes_list=vector("list",length=10)
for(c in 1:10){
  #SNPs
  geno=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
  fes=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_male_flowering_d6_x_%s_founderprobs.rds',c,env))
  lowrep=readRDS(sprintf('genotypes/probabilities/geno_probs/dropped/bg%.0f_low_rep_markers.rds',c))
  betas = fes[,grep('beta',colnames(fes))]
  rownames(betas)=fes$X_ID
  betas= betas[complete.cases(betas),]
  betas[,-1] = betas[,-1]+betas[,1]
  names(lowrep)=paste0('beta.',seq(1,16))
  for(i in 1:16){
    b=paste0('beta.',i)
    lowm=lowrep[[b]]
    lowm=lowm[lowm %in% rownames(betas)]
    if(!is.null(lowm)){
      betas[lowm,b]=0
    }
    betas[betas[,b] < 0.05 & betas[,b] > -0.05,b] = 0
  }


  #rownames(geno)=geno$ind
  #geno=geno[inds,fes$X_ID]
  geno_list[[c]]=geno
  fes_list[[c]]=betas
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
  betas=fes_list[[c]]
  nmarkers=nrow(betas)
  markers=rownames(betas)
  t=matrix(unlist(lapply(geno,function(x) x[ind,markers])),nrow=nmarkers)
  colnames(t)=names(geno)
  rownames(t)=markers
  pg=sapply(markers,function(m) sum(t[m,] * betas[m,]))
  #pg=sapply(seq(1,length(t)), function(x) betas[names(t)[x],b[x]])
  #pg=sapply(seq(1,nmarkers), function(i) betas[rownames(betas)[i],paste0('beta.',(geno[ind,)+1)])
  pg=sum(pg)
  names(pg)=ind
  return(pg)
}
#print(system.time({
#  results=chr_pgs(10,'EB.09S.H.00002')
#}))


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
print(system.time({
  results=calc_pgs('EB.09S.H.00002')
}))

n_inds=inds

print(system.time({
results=mclapply(n_inds,calc_pgs,mc.cores=cores)
}))

all_pgs=data.frame(inds=inds,f_pgs=unlist(results),stringsAsFactors=F)
fwrite(all_pgs,'GridLMM/pgs/founder_dta_polygenic_scores.txt',quote=F,sep='\t',row.names=F)
