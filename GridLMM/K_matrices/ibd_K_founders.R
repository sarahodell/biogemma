#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')

X_list=readRDS(sprintf('../../geno_probs/bg%s_filtered_genotype_probs.rds',c))
gmap=fread(sprintf('../../qtl2_startfiles/Biogemma_gmap_c%s.csv',c),data.table=F,stringsAsFactors=F)
pmap=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F,stringsAsFactors=F)

map=merge(gmap,pmap[,c('marker','pos')],by="marker")
names(map)=c('marker','chr','cM','pos')

markers=c(dimnames(X_list[[1]])[[2]])

submap=map[map$marker %in% markers,]
submap=submap[order(submap$pos),]
rownames(submap)=seq(1,dim(submap)[1])
submap$cdist=c(0,diff(as.matrix(submap$cM)))
markers_ord=c(submap$marker)

X_list_ordered=lapply(X_list,function(i) i[,markers_ord])

X_list_full = lapply(X_list,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.95,1,ifelse(x[,i]<=0.05,0,round(x[,i],2)))))
for(i in 1:length(X_list_full)){dimnames(X_list_full[[i]])[[2]]=dimnames(X_list_ordered[[i]])[[2]]}

shared=array(0,dim=c(344,344))

for(i in 1:length(X_list_full)){
  print(i)
  t=X_list_full[[i]]
  markers=dimnames(t)[[2]]
  for(x in 1:dim(t)[2]){
    m=t[,x]
    m_array=matrix(0,nrow=344,ncol=344)
    has=which(m==1,m)
    if(length(has)!=0){
      m_array[,has]=m
      m_array=m_array * submap[submap$marker==markers[x],]$cdist
      shared=shared+m_array
    }
  }
}
for(i in 1:dim(shared)[1]){shared[i,i]=max(map$cM)}
dimnames(shared)=list(dimnames(X_list[[1]])[[1]],dimnames(X_list_ordered[[1]])[[1]])


fwrite(shared,sprintf('chr%s_founder_ibd_K_matrix.txt',c),quote=F,row.names=F,col.names=T)

