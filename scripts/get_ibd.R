#!/usr/bin/env Rscript

library("qtl2")
library("data.table")
library("ggplot2")
library("reshape2")
library("tibble")
library('igraph')
library('abind')

args=commandArgs(trailingOnly=T)
c=as.character(args[1])

print(c)
bg=read_cross2(sprintf('../Biogemma_071918/start_files/Biogemma_c%s.json',c))
ibd=find_ibd_segments(bg$founder_geno,bg$pmap,min_lod=15,error_prob=0.002,cores=4)
fwrite(ibd,sprintf('bg%s_pairwise_ibdsegments_121018.txt',c))

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
founderfile=fread('Founder_colorcodes.txt',header=F)
names(founderfile)=c('line','hex')
code=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
founderfile$code=code
i1 <- match(ibd$strain1,founderfile$code)
i2 <- match(ibd$strain2,founderfile$code)
ibd$line1=founderfile$line[i1]
ibd$line2=founderfile$line[i2]

ibd$pair=do.call(paste,c(ibd[c('line1','line2')],sep="-"))

ibd=ibd[order(ibd$left_pos),]
rownames(ibd)=seq(1,nrow(ibd))

ibd_segments=c()
start=min(ibd$left_pos) #start with first segment in chromosome
ibd_graph=array(0,dim=c(16,16,1))
while(start<max(ibd$right_pos)){
  remaining=ibd[(ibd$left_pos>start | ibd$right_pos>=start),] #grab segments that are to the right of the start
  remaining$left_dist=remaining$left_pos - start 
  remaining$right_dist=remaining$right_pos - start
  lowest_left=min(remaining$left_dist[remaining$left_dist>0])
  lowest_right=min(remaining$right_dist[remaining$right_dist>0])
  closest=min(lowest_left,lowest_right) #find the next start or end of segment that is closest to current start
   if (lowest_left==closest){
    next_break=remaining[remaining$left_dist==closest,]$left_pos[1]
  }
  if (lowest_right==closest){
    next_break=remaining[remaining$right_dist==closest,]$right_pos[1]
  }
  within=remaining[(remaining$right_pos>=next_break) & (remaining$left_pos<=start),] #identify ibd segments that are between start and next_break

  if(dim(within)[1]>=1){
    adj=matrix(0,nrow=16,ncol=16,dimnames=list(founders,founders))
    for (i in 1:nrow(within)){
      adj[c(within[i,'line1']),c(within[i,'line2'])]=1
      adj[c(within[i,'line2']),c(within[i,'line1'])]=1
    }
    ibd_graph=abind(ibd_graph,adj)
    graph=graph_from_adjacency_matrix(adj,mode=c('undirected'))
    blocks=unname(components(graph)$membership)
    n_grps=components(graph)$no
  }
  else{
    blocks=seq(1,16)
    n_grps=16
  }
  ibd_segments=rbind(ibd_segments,c(10,start,next_break,blocks,n_grps))
  start=next_break
}
ibd_segments=as.data.frame(ibd_segments)
names(ibd_segments)=c('chrom','start','end',founders,'n_haps')
ibd_segments$start=ibd_segments$start*1e6
ibd_segments$end=ibd_segments$end*1e6
saveRDS(ibd_graph,sprintf('bg%s_ibd_graph.rds',c))
fwrite(ibd_segments,sprintf('bg%s_ibd_blocks.txt',c),row.names=F,quote=F,sep='\t')

pr=readRDS(sprintf('../Biogemma_082318/600K_DHgenoprobs/bg%s_0823_genoprobs.rds',c))
pmap=fread(sprintf('../Biogemma_082318/Biogemma_pmap_c%s.csv',c),data.table=F)
pmap$pos=pmap$pos*1e6
dimnames(pr[[1]])[[2]]=founders
samples=unlist(dimnames(pr[[1]])[1])

final_haplo=list(chr=c)

groups=sort(unique(ibd_segments$n_haps))
for(h in groups){
  n_hap=h
  n_ind=dim(pr[[1]])[1]
  # number of IBD Segments with 15 unique haplotypes
  test=ibd_segments[ibd_segments$n_haps==n_hap,]
  rownames(test)=seq(1,nrow(test))
  # for each individual
  for(n in seq(1,n_ind)){
    hprobs=c()
    markers=c()
    #go through each haplotype block with n_hap haplotype groups
    for(i in seq(1,nrow(test))){
      start=test$start[i]
      end=test$end[i]
      # grab the SNPs within these segments
      within=pmap[(pmap$pos < end) & (pmap$pos >= start),]
      if(dim(within)[1]!=0){
        # get the founder probabilities of these SNPs
        prob=pr[[1]][n,,within$marker]
        line=test[i,]
        markers=c(markers,within$marker)
        # for each haplotype group grab the founders in that haplotype group
        # and sum the probability
        #if(length(dim(prob))==2){
        hprob=sapply(seq(1,n_hap), function(x) t(prob) %*% as.vector(line[,founders]==x))
        #else{
        #  hap_prob=sapply(seq(1,n_hap), function(x) t(prob[n,,]) %*% as.vector(line[,founders]==x))
        #}
        hprobs=rbind(hprobs,hprob)
      }
    }
    a=dim(hprobs)[1]
    b=dim(hprobs)[2]
    hprob_array=array(hprobs,dim=c(1,a,b))
    if(n==1){
      haplo_probs=array(hprob_array,dim=c(n_ind,a,b))
    }
    else{
      haplo_probs[n,,]=hprob_array
    }
  }
  dimnames(haplo_probs)=list(samples,markers,seq(1,h))
  final_haplo[[h]]=haplo_probs
  
}

# for each IBD segment
saveRDS(final_haplo,sprintf('bg%s_haplotype_probs.rds',c))
