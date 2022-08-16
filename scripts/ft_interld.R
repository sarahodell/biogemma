#!/usr/bin/env Rscript

library('data.table')
library('reshape2')
library('dplyr')

ftbed=fread('selection/FT_gene_list_AGPv4.bed',data.table=F,header=F)
names(ftbed)=c('chrom','start','end','geneID')

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
#pos=unique(c(paste0(interld$CHR_A,'_',interld$BP_A),paste0(interld$CHR_B,'_',interld$BP_B)))

ft_counts=c()
overlaps=c()

for(c in 1:10){
  ftsub=ftbed[ftbed$chrom==c,]
  rownames(ftsub)=seq(1,nrow(ftsub))
  subpos=pos[pos$chr==c,]

  env1=subpos
  env1=as.data.table(subpos)
  env2=as.data.table(ftsub)
  setkey(env2,start,end)
  comparison=foverlaps(env1,env2,by.x=c('start','end'),by.y=c('start','end'),nomatch=NULL)
  count=c(c,length(unique(comparison$geneID)))
  ft_counts=rbind(ft_counts,count)
  overlaps=rbind(overlaps,comparison)
}

ft_counts=as.data.frame(ft_counts,stringsAsFactors=F)
rownames(ft_counts)=seq(1,nrow(ft_counts))
names(ft_counts)=c('chrom','overlapping_gene_count')
sum(ft_counts$overlapping_gene_count)
#[1] 593

fwrite(ft_counts,'selection/interchrom_ld_ft_gene_counts.txt',quote=F,sep='\t',row.names=F)
