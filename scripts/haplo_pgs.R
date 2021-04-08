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
#s_pgs=rep(0,length(inds))
#f_pgs=rep(0,length(inds))
#h_pgs=rep(0,length(inds))

geno_list=vector("list",length=10)
hes_list=vector("list",length=10)
for(c in 1:10){
  base=base_list[c]
  count=1
  hapnames=c()
  for(h in base:16){
    hapfile=sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%.0f_filtered_haplogroup%.0f_probs.rds',c,h)
    hesfile=sprintf('GridLMM/GridLMM_haplotypes/models/Biogemma_chr%.0f_haplogrp%.0f_male_flowering_d6_x_%s.rds',c,h,env)
    if(file.exists(hapfile) & file.exists(hesfile)){
      geno=readRDS(hapfile)
      hes=readRDS(hesfile)
      #lowrep=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/'))
      lowrep=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/low_rep/bg%.0f_haplogroup%.0f_low_rep_markers.rds',c,h))
      betas = hes[,grep('beta',colnames(hes))]
      rownames(betas)=hes$X_ID
      betas= betas[complete.cases(betas),]
      names(lowrep)=paste0('beta.',seq(1,h))
      betas[,-1] = betas[,-1]+betas[,1]
      for(i in 1:h){
        b=paste0('beta.',i)
        lowm=lowrep[[b]]
        lowm=lowm[lowm %in% rownames(betas)]
        if(!is.null(lowm)){
          betas[lowm,b]=0
        }
        betas[betas[,b] < 0.05 & betas[,b] > -0.05,b] = 0
      }
      #pg=mapply(function(i,j) betas[i,j], grab_row, grab_col)
      #rownames(geno)=geno$ind
      #geno=geno[inds,fes$X_ID]
      geno_list[[c]][[count]]=geno
      hes_list[[c]][[count]]=betas
      hapnames=c(hapnames,as.character(h))
      count=count+1
    }
  }
  names(hes_list[[c]])=hapnames
  names(geno_list[[c]])=hapnames
}

hap_pgs<-function(h,c,ind){
  geno=geno_list[[c]][h][[1]]
  bts=as.data.frame(hes_list[[c]][h],stringsAsFactors=F)
  nmarkers=nrow(bts)
  markers=rownames(bts)
  t=matrix(unlist(lapply(geno,function(x) x[ind,markers])),nrow=nmarkers)
  colnames(t)=names(geno)
  rownames(t)=markers
  pg=sapply(markers,function(m) sum(t[m,] * bts[m,]))
  pg=sum(pg)
}

chr_pgs<-function(c,ind){
  base=base_list[c]
  pg=0
  for(h in base:16){
    h=as.character(h)
    if(h %in% names(hes_list[[c]])){
      pg=pg+hap_pgs(h,c,ind)
    }
  }
  #pg=sapply(seq(1,length(t)), function(x) betas[names(t)[x],b[x]])
  #pg=sapply(seq(1,nmarkers), function(i) betas[rownames(betas)[i],paste0('beta.',(geno[ind,)+1)])
  names(pg)=ind
  return(pg)
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
print(system.time({
  results=calc_pgs('EB.09S.H.00002')
}))

n_inds=inds

print(system.time({
results=mclapply(n_inds,calc_pgs,mc.cores=cores)
}))

all_pgs=data.frame(inds=inds,h_pgs=unlist(results),stringsAsFactors=F)
fwrite(all_pgs,'GridLMM/pgs/haplotype_dta_polygenic_scores.txt',quote=F,sep='\t',row.names=F)
