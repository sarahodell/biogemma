#!/usr/bin/env Rscript

library('data.table')

full_length=0

full_k=array(0,dim=c(344,344))
for(c in 1:10){
  #hk=array(0,dim=c(344,344))
  #for(h in 2:16){
    #dh=fread(sprintf('chr%.0f_hapgrp%.0f_ibd_K_matrix.txt',c,h),data.table=F)
  #  names=names(dh)
  #  dh=as.matrix(dh)
  #  chrlen=dh[1,1]
  #  #for(i in 1:dim(dh)[1]){dh[i,i]=0}
  #  hk=hk+dh
  #}
  dh=fread(sprintf('chr%.0f_founder_ibd_K_matrix.txt',c),data.table=F,stringsAsFactors=F)
  chrlen=dh[1,1]
  names=names(dh)
  dh=as.matrix(dh)
  full_length=full_length+chrlen
  full_k=full_k + dh
}
print(full_length)
#for(i in 1:dim(full_k)[1]){full_k[i,i]=full_length}

full_k_prop=full_k/full_length
#fwrite(full_k_prop,'full_founder_ibd_kinship.txt',quote=F,col.names=T)
dimnames(full_k_prop)=list(names,names)
fwrite(full_k_prop,'full_founder_ibd_kinship.txt',quote=F,col.names=T)