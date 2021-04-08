#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])
h=as.numeric(args[[2]])

library("data.table")
library("readr")
library("ggplot2")
library("reshape2")
library("tibble")
library('igraph')
library('abind')

#founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
pmapfile=sprintf("genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv",c)
ibdfile=sprintf('ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',c)

# Read in IBD segments from get_ibd.R

pmap=fread(pmapfile,data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

hapfile=sprintf("genotypes/probabilities/haplotype_probs/RefinedIBD_600K/raw/bg%s_refined_ibd_haplogroup%.0f_probs.rds",c,h)
ibd_segments=fread(ibdfile,data.table=F)
if(file.exists(hapfile)){
  haplo_probs=readRDS(hapfile)
  size=dim(haplo_probs[[1]])[2]
  f_sums=lapply(haplo_probs,colSums)
  low_rep=lapply(f_sums,function(x) which(x<.99,arr.ind=T))

  dropped=vector("list",length=h)
  for(i in 1:h){
    drop=unlist(unname(low_rep[[i]]))
    if(length(drop)!=0){
      dropped[[i]]=dimnames(haplo_probs[[i]])[[2]][drop]
      #print(i)
      #print(length(drop))
      #haplo_probs[[i]]=haplo_probs[[i]][,drop]=NA
    }
  }


  # drop sites with less than 5 lines with prob greater than .8
  #Remove sites with low haplotype representation
  # drop sites with summed founder prob of less than 1
  f_sums2=lapply(haplo_probs, function(x) which(colSums(x>=0.8)<5,arr.ind=T))

  for(i in 1:h){
    drop=unlist(unname(f_sums2[[i]]))
    if(length(drop)!=0){
      dropped[[i]]=unique(c(dropped[[i]],dimnames(haplo_probs[[i]])[[2]][drop]))
    }
  }

  saveRDS(dropped,sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/low_rep/bg%s_haplogroup%.0f_low_rep_markers.rds',c,h))

  #regrp=c()
  #for(i in 1:h){regrp[[i]]=haplo_probs[[1]][,,i]}
  #r2 cutoff of 0.95
  cutoff=0.95

  # Filter Correlated Markers in Haplotype Probability
  m_names=names(haplo_probs[[1]][1,])
  dropped=list()
  count=1
  dropped[[count]]=list(marker=c(m_names[1]),linked=c())

  keep=c(1)
  start=c()
  for(x in 1:h){start=c(start,as.vector(haplo_probs[[x]][,1]))}
  size=dim(haplo_probs[[1]])[2]
  for(i in 1:size){
    ind_max = c()
    for(x in 1:h){ind_max=c(ind_max,as.vector(haplo_probs[[x]][,i]))}
    if((cor(start,ind_max,use="complete.obs")**2)<cutoff){
      keep=c(keep,i)
      start=ind_max
      start_ind=i
      count=count+1
      dropped[[count]]=list(marker=c(m_names[start_ind]),linked=c())
    }
    else{
      dropped[[count]]$linked=c(dropped[[count]]$linked,m_names[i])
    }
  }

  saveRDS(dropped,sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/dropped/bg%s_haplogroup%s_dropped_markers.rds',c,h))
  filtered_hap=c()
  for(k in 1:h){filtered_hap[[k]]=haplo_probs[[k]][,keep]}
  if(is.null(dim(filtered_hap[[1]]))){
    for(k in 1:h){
      dimk=length(filtered_hap[[k]])
      filtered_hap[[k]]=array(filtered_hap[[k]],dim=c(dimk,1),dimnames=list(names(filtered_hap[[k]]),m_names[keep]))
      }
  }
  saveRDS(filtered_hap,sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%s_probs.rds',c,h))
  print("Finished filtering correlated markers")
  sprintf("Keeping %s markers in haplogroup %s on chromosome %s with correlation filter of %.2f",length(keep),h,c,cutoff)
  print(sprintf("Finished haplotype %s",h))
}else{
  print(sprintf("No markers in %s haplogroup %.0f",c,h))
}
