#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')

fgeno=fread(sprintf('../Biogemma_foundergenos/Founder_genos_chr%s_121718.csv',c),data.table=F)
pmap=fread(sprintf('qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

ibd=fread(sprintf('ibd_segments/bg%s_ibd_segments_adjusted.txt',c),data.table=F)
old_ibd=fread('sprintf('ibd_segments/bg%s_wgs_ibdsegments_010319.txt',c),data.table=F)

print(dim(ibd))
print(dim(old_ibd))

mismatches=c()
mismcatch2=c()

for(i in seq(1,dim(ibd)[1])){
      line=ibd[i,]
      old_line=old_ibd[i,]
      window=fgeno[c(line$strain1,line$strain2),line$left_index:line$right_index]
      mm=sum(apply(window,MARGIN=2,FUN=function(x) length(unique(x)))>1)
      if(mm > line$n_mismatch){
      	    mismatches=c(mismatches,i)
      }
      if(mm > old_line$n_mismatch){
      	    mismcatch2=c(mismatch2,i)

}

print(dim(ibd)[1])
print(length(mismatches))
print(length(mismatch2))