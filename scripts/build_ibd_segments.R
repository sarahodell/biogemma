#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library("data.table")
library("readr")
library("ggplot2")
library("reshape2")
library("tibble")
library('igraph')
library('abind')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
pmapfile=sprintf("genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv",c)
ibd=fread(sprintf('ibd_segments/refinedibd/600K/Biogemma_600K_Founders_RefinedIBD_chr%s.ibd',c),data.table = F)

names(ibd)=c('ID_1','Hap_ID1','ID_2','Hap_ID2','Chromosome','left_pos','right_pos','LOD','cM_length')
ibd=ibd[ibd$ID_1!="MBS847" & ibd$ID_2 !="MBS847",]
ibd=ibd[order(ibd$left_pos),]
rownames(ibd)=seq(1,dim(ibd)[1])

print("Finished formating IBD file")

pmap=fread(pmapfile,data.table=F)
#hap_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

ibd_segments=c()
start=min(ibd$left_pos) #start with first segment in chromosome
ibd_graph=array(0,dim=c(16,16,1))
ibd_graph[,,1]=diag(16)
dimnames(ibd_graph)=list(founders,founders,"blank")
if(start>min(pmap$pos)){
  ibd_segments=rbind(ibd_segments,c(c,min(pmap$pos),start,seq(1,16),16))
}

while(start<max(ibd$right_pos)){
  remaining=ibd[(ibd$left_pos>start | ibd$right_pos>=start),] #grab segments that are to the right of the start
  remaining$left_dist=remaining$left_pos - start
  remaining$right_dist=remaining$right_pos - start
  if(length(remaining$left_dist[remaining$left_dist>0]!=0)){
    lowest_left=min(remaining$left_dist[remaining$left_dist>0])
  }
  else{
    lowest_left=Inf
  }
  if(length(remaining$right_dist[remaining$right_dist>0])!=0){
    lowest_right=min(remaining$right_dist[remaining$right_dist>0])
  }
  else{
    lowest_right=Inf
  }
  closest=min(lowest_left,lowest_right) #find the next start or end of segment that is closest to current start
  if (lowest_left==closest){
    next_break=remaining[remaining$left_dist==closest,]$left_pos[1]
  }
  if (lowest_right==closest){
    next_break=remaining[remaining$right_dist==closest,]$right_pos[1]
  }
  within=remaining[(remaining$right_pos>=next_break) & (remaining$left_pos<=start),] #identify ibd segments that are between start and next_break

  if(dim(within)[1]>=1){
    rownames(within)=seq(1,dim(within)[1])
    #adj=matrix(0,nrow=16,ncol=16,dimnames=list(founders,founders))
    adj=diag(16)
    dimnames(adj)=list(founders,founders)
    for (i in 1:nrow(within)){
      adj[c(within[i,'ID_1']),c(within[i,'ID_2'])]=1
      adj[c(within[i,'ID_2']),c(within[i,'ID_1'])]=1
    }
    ibd_graph=abind(ibd_graph,adj)
    graph=graph_from_adjacency_matrix(adj,mode=c('undirected'))
    blocks=unname(components(graph)$membership)
    n_grps=components(graph)$no
  }
  else{
    adj=diag(16)
    dimnames(adj)=list(founders,founders)
    ibd_graph=abind(ibd_graph,adj)
    blocks=seq(1,16)
    n_grps=16
  }
  ibd_segments=rbind(ibd_segments,c(as.integer(c),as.integer(start),as.integer(next_break),as.integer(blocks),as.integer(n_grps)))
  start=next_break
}


if(start<max(pmap$pos)){
  ibd_segments=rbind(ibd_segments,c(c,start,max(pmap$pos),seq(1,16),16))
}

ibd_segments=as.data.frame(ibd_segments,stringsAsFactors=F)
names(ibd_segments)=c('chrom','start','end', founders,'n_haps')
print("Finished making IBD block file")

saveRDS(ibd_graph,sprintf('ibd_segments/refinedibd/600K/bg%s_ibd_graph.rds',c))
fwrite(ibd_segments,sprintf('ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',c),row.names=F,quote=F,sep='\t')

groups=sort(unique(as.integer(ibd_segments$n_haps)))

line=data.table(chr=as.numeric(c),min_hap=min(groups))
fwrite(line,file="genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt",sep='\t',append=T,col.names=F,quote=F,row.names=F)
