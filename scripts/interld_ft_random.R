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

interld=fread('stats/ld_decay/Biogemma_DHLines_rsquared_all_chroms_r2_0.9.ld',data.table=F)
pos=interld[,c('CHR_A','BP_A')]
pos2=interld[,c('CHR_B','BP_B')]
names(pos2)=c('CHR_A','BP_A')
pos=rbind(pos,pos2)
pos=pos[,c('CHR_A','BP_A','BP_A')]
#pos$name=paste0(pos$CHR_A,'_',pos$BP_A)
pos=pos[!duplicated(pos), ]
names(pos)=c('chr','start','end')
rownames(pos)=seq(1,nrow(pos))

total_segs=0
seg_size=0

for(c in 1:10){
  randsub=randbed[randbed$chrom==c,]
  rownames(randsub)=seq(1,nrow(randsub))
  subpos=pos[pos$chr==c,]
  total_segs=total_segs+nrow(subpos)
  #names(chi_segs)=c('chrom','start','end','markers','hapgrp')
  #seg_size=seg_size+sum(chi_segs$end-chi_segs$start)
  env1=subpos
  env1=as.data.table(env1)
  env2=as.data.table(randsub)
  setkey(env2,start,end)
  comparison=foverlaps(env1,env2,by.x=c('start','end'),by.y=c('start','end'),nomatch=NULL)
  count=c(c,length(unique(comparison$geneID)))
  rep_counts=rbind(rep_counts,count)
}

rep_counts=as.data.frame(rep_counts,stringsAsFactors=F)
rownames(rep_counts)=seq(1,nrow(rep_counts))
names(rep_counts)=c('chrom','overlapping_gene_count')

fwrite(rep_counts,sprintf('selection/reps/interchrom_ld_random_gene_counts_rep%s.txt',rep),quote=F,sep='\t',row.names=F)
