#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
rep=as.character(args[[1]])

library('data.table')
library('tidyverse')

ftbed=fread('selection/FT_gene_list_AGPv4.bed',data.table=F,header=F)
all_bed=fread('selection/Zea_mays_AGPv4_full_gene_list.bed',data.table=F,header=F)
size=nrow(ftbed)
sample=sample(seq(1,nrow(all_bed)),size)

randbed=all_bed[sample,]
rownames(randbed)=seq(1,size)
names(randbed)=c('chrom','start','end','geneID')
rep_counts=c()

baselist=c(8,7,7,8,6,7,7,8,7,7)
total_segs=0
seg_size=0

for(c in 1:10){
  base=baselist[c]
  randsub=randbed[randbed$chrom==c,]
  rownames(randsub)=seq(1,nrow(randsub))
  chi_segs=fread(sprintf('selection/haplotype_probs/bg%.0f_chi_peak_seqments.bed',c),data.table=F,header=F)
  total_segs=total_segs+nrow(chi_segs)
  names(chi_segs)=c('chrom','start','end','markers','hapgrp')
  seg_size=seg_size+sum(chi_segs$end-chi_segs$start)
  env1=chi_segs
  env1=as.data.table(env1)
  env2=as.data.table(randsub)
  setkey(env2,start,end)
  comparison=foverlaps(env1,env2,by.x=c('start','end'),by.y=c('start','end'),nomatch=NULL)
  within=nrow(comparison[(comparison$start>=comparison$i.start & comparison$end <= comparison$i.end),])
  hapcounts=comparison %>% group_by(hapgrp) %>% count
  all_hapcounts=c()
  #hapgrps=paste0('hapgrp',seq(6,16))
  base=baselist[c]
  for(h in 6:16){
    hap=paste0('hapgrp',h)
    if(base>h){
      all_hapcounts=c(all_hapcounts,NA)
    }
    else if(hap %in% hapcounts$hapgrp){
      all_hapcounts=c(all_hapcounts,hapcounts[hapcounts$hapgrp==hap,]$n)
    }
    else{
      all_hapcounts=c(all_hapcounts,0)
    }
  }
  count=c(c,nrow(comparison),within,all_hapcounts)
  rep_counts=rbind(rep_counts,count)
}

rep_counts=as.data.frame(rep_counts,stringsAsFactors=F)
rownames(rep_counts)=seq(1,nrow(rep_counts))
names(rep_counts)=c('chrom','overlapping_gene_count','within_gene_count',paste0('hapgrp',seq(6,16)))

fwrite(rep_counts,sprintf('selection/reps/haplotype_random_gene_counts_rep%s.txt',rep),quote=F,sep='\t',row.names=F)
