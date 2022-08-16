#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('tidyverse')

baselist=c(8,7,7,8,6,7,7,8,7,7)
min=baselist[as.numeric(c)]

count=1
len=length(min:16)
markers=vector("list",length=len)

for(i in min:16){
  #mod=readRDS(sprintf('raw/bg%s_refined_ibd_haplogroup%.0f_probs.rds',c,i))
  mod=readRDS(sprintf('bg%s_filtered_haplogroup%.0f_probs.rds',c,i))
  markers[[count]]=dimnames(mod[[1]])[[2]]
  count=count+1
}

all_markers=c()
for(i in 1:len){
  all_markers=c(all_markers,markers[[i]])
}
length(all_markers)

length(unique(all_markers))

if(length(all_markers)>length(unique(all_markers))){
  print(sprintf("Duplicate Markers on chromosome %s",c))
  for(i in 1:(len-1)){
    for(j in (i+1):len){
      both=c(markers[[i]],markers[[j]])
      if(length(both)>length(unique(both))){
        print(sprintf("Duplicate Markers in %s and %s",seq(min,16)[i],seq(min,16)[j]))
      }
    }
  }
}

low_rep=c()
for(i in min:16){
  #mod=readRDS(sprintf('raw/bg%s_refined_ibd_haplogroup%.0f_probs.rds',c,i))
  mod=readRDS(sprintf('low_rep/bg%s_haplogroup%.0f_low_rep_markers.rds',c,i))
  for(j in 1:min){
    low_rep=unique(c(low_rep,mod[[j]]))
  }
}
print("Number of filtered markers that have low representation")
print(sum(all_markers %in% low_rep))
