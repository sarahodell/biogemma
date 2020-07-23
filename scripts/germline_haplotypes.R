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

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#founders2=c("A632","B73","CO255","FV252","OH43", "A654","FV2","C103","EP1","D105","W117","B96","DK63","F492","ND245","VA85")

ibd=fread(sprintf('ibd_segments/refinedibd/600K/Biogemma_600K_Founders_RefinedIBD_chr%s.ibd',c),data.table = F)
#ibd=read_table2(sprintf('ibd_segments/germline/600K/Biogemma_Founders_germline_IBD_chr%s.match',c),col_names = F)
#ibd=as.data.frame(ibd,stringsAsFactors=F)
#for(i in seq(1,dim(ibd)[1])){
#      if(!(ibd[i,]$X1 %in% founders)){
#      	f2=which(founders2==ibd[i,]$X1)
#        ibd[i,]$X1=founders[f2]
#
#	f2=which(founders2==ibd[i,]$X3)
#  	ibd[i,]$X3=founders[f2]
#  }
#}

#ibd=ibd[,c(1,3,5:15)]

#names(ibd)=c('ID_1','Family_ID_1','ID_2','Family_ID2','Chromosome','left_pos','right_pos','left_marker','right_marker','n_markers','Genetic_length','Genetic_length_units','n_mismatch','ID1_Is_Homozygous','ID2_Is_Homozygous')
names(ibd)=c('ID_1','Hap_ID1','ID_2','Hap_ID2','Chromosome','left_pos','right_pos','LOD','cM_length')
ibd=ibd[ibd$ID_1!="MBS847" & ibd$ID_2 !="MBS847",]

#ibd=ibd[ibd$ID_1!="MBS847" & ibd$ID_2!="MBS847",]
ibd=ibd[order(ibd$left_pos),]
rownames(ibd)=seq(1,dim(ibd)[1])

print("Finished formating IBD file")


genofile=sprintf("genotypes/probabilities/geno_probs/raw/bg%s_genoprobs_010319.rds",c)
pmapfile=sprintf("genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv",c)
outfile=sprintf("genotypes/probabilities/haplotype_probs/bg%s_refined_ibd_haplotype_probs.rds",c)

# Read in IBD segments from get_ibd.R

pmap=fread(pmapfile,data.table=F)
#founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
hap_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

ibd_segments=c()
start=min(ibd$left_pos) #start with first segment in chromosome
ibd_graph=array(0,dim=c(16,16,1))
ibd_graph[,,1]=diag(16)
dimnames(ibd_graph)=list(hap_founders,hap_founders,"blank")
if(start>min(pmap$pos)){
  ibd_segments=rbind(ibd_segments,c(c,min(pmap$pos),start,seq(1,16),16))
}

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
    rownames(within)=seq(1,dim(within)[1])
    #adj=matrix(0,nrow=16,ncol=16,dimnames=list(founders,founders))
    adj=diag(16)
    dimnames(adj)=list(hap_founders,hap_founders)
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
    dimnames(adj)=list(hap_founders,hap_founders)
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
names(ibd_segments)=c('chrom','start','end',hap_founders,'n_haps')
#dimnames(ibd_graph)[[3]]=c("blank",ibd_segments$start)
print("Finished making IBD block file")

saveRDS(ibd_graph,sprintf('ibd_segments/refinedibd/600K/bg%s_ibd_graph.rds',c))
fwrite(ibd_segments,sprintf('ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',c),row.names=F,quote=F,sep='\t')

pr=readRDS(genofile)
samples=unlist(dimnames(pr[[1]])[1])

groups=sort(unique(as.integer(ibd_segments$n_haps)))

line=data.table(chr=as.numeric(c),min_hap=min(groups))
fwrite(line,file="genotypes/probabilities/haplotype_probs/min_haps.txt",sep='\t',append=T,col.names=F,quote=F,row.names=F)

final_haplo=vector("list",length=length(groups))
names(final_haplo)=groups
for(h in groups){
  print(sprintf("Starting haplotype %s",h))
  n_hap=h
  n_ind=dim(pr[[1]])[1]
  # number of IBD Segments with 15 unique haplotypes
  test=ibd_segments[ibd_segments$n_haps==n_hap,]
  test=test[order(test$start),]
  rownames(test)=seq(1,nrow(test))
  haplo_probs=NULL
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
      if(dim(within)[1]!=0 & is.null(dim(within))==F){
        # get the founder probabilities of these SNPs
        prob=pr[[1]][n,,within$marker]
        line=test[i,]
        markers=c(markers,within$marker)
        # for each haplotype group grab the founders in that haplotype group
        # and sum the probability
        hprob=sapply(seq(1,n_hap), function(x) t(prob) %*% as.vector(line[,hap_founders]==x))
        hprobs=rbind(hprobs,hprob)
      }
    }
    a=dim(hprobs)[1]
    b=dim(hprobs)[2]
    hprob_array=array(hprobs,dim=c(1,a,b))
    if(is.null(haplo_probs)){
      haplo_probs=hprob_array
    }
    else{
      haplo_probs=abind(haplo_probs,hprob_array,along=1)
    }
  }
  dimnames(haplo_probs)=list(samples,markers,seq(1,h))
# final_haplo[[h]]=haplo_probs
  regrp=c()
  for(i in 1:h){regrp[[i]]=haplo_probs[,,i]}
  cutoff=0.95

  # Filter Correlated Markers in Haplotype Probability
  m_names=names(regrp[[1]][1,])
  dropped=list()
  count=1
  dropped[[count]]=list(marker=c(m_names[1]),linked=c())

  keep=c(1)
  start=c()
  for(x in 1:h){start=c(start,as.vector(regrp[[x]][,1]))}
  for(i in 2:dim(regrp[[1]])[2]){
    ind_max = c()
    for(x in 1:h){ind_max=c(ind_max,as.vector(regrp[[x]][,i]))}
    if(cor(start,ind_max)<cutoff){
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

  saveRDS(dropped,sprintf('genotypes/probabilities/haplotype_probs/bg%s_haplogroup%s_dropped_markers.rds',c,h))
  filtered_hap=c()
  for(k in 1:h){filtered_hap[[k]]=regrp[[k]][,keep]}
  if(is.null(dim(filtered_hap[[1]]))){
    for(k in 1:h){
      dimk=length(filtered_hap[[k]])
      filtered_hap[[k]]=array(filtered_hap[[k]],dim=c(dimk,1),dimnames=list(names(filtered_hap[[k]]),m_names[keep]))
      }
  }
  saveRDS(filtered_hap,sprintf('genotypes/probabilities/haplotype_probs/bg%s_filtered_haplogroup%s_probs_2.rds',c,h))
  print("Finished filtering correlated markers")
  sprintf("Keeping %s markers in haplogroup %s on chromosome %s with correlation filter of %.2f",length(keep),h,c,cutoff)
  #final_haplo[[as.character(h)]]=haplo_probs
  print(sprintf("Finished haplotype %s",h))
}


#saveRDS(final_haplo,outfile)
