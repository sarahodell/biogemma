#!/usr/bin/env Rscript


library('data.table')
library('tidyverse')

ftbed=fread('selection/FT_gene_list_AGPv4.bed',data.table=F,header=F)
names(ftbed)=c('chrom','start','end','geneID')

ft_counts=c()

baselist=c(8,7,7,8,6,7,7,8,7,7)

for(c in 1:10){
  base=baselist[c]
  ftsub=ftbed[ftbed$chrom==c,]
  rownames(ftsub)=seq(1,nrow(ftsub))
  chi_segs=fread(sprintf('selection/haplotype_probs/bg%.0f_chi_peak_seqments.bed',c),data.table=F,header=F)
  names(chi_segs)=c('chrom','start','end','markers','hapgrp')
  env1=chi_segs
  env1=as.data.table(env1)
  env2=as.data.table(ftsub)
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
  ft_counts=rbind(ft_counts,count)
}

ft_counts=as.data.frame(ft_counts,stringsAsFactors=F)
rownames(ft_counts)=seq(1,nrow(ft_counts))
names(ft_counts)=c('chrom','overlapping_gene_count','within_gene_count',paste0('hapgrp',seq(6,16)))

fwrite(ft_counts,'selection/haplotype_ft_gene_counts.txt',quote=F,sep='\t',row.names=F)
