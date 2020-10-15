#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')
library('ggplot2')
library('tidyverse')

total_tests=125764
bonf=-log10(0.05 / total_tests)
hap_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=T)
gmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_gmap_c%s.csv',c),data.table=T)

ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',c),data.table=T)

get_weighted_null=function(i,comparison){
  weights=melt(comparison[i,c('marker',..hap_founders)],'marker') %>% count(value)
  null=rep(1/16,h)
  weighted_null=weights$n*null
  return(weighted_null)
}

haplo_counts=list()
count=0
all_p_chi=c()
baselist=c(8,7,7,8,6,7,7,8,7,7)
for(h in baselist[as.numeric(c)]:16){
  test=ibd[ibd$n_haps==h,]
  rownames(test)=seq(1,dim(test)[1])

  env1=pmap
  env1$pos_end=env1$pos
  env1=as.data.table(env1)
  env2=test
  env2$end=env2$end-1
  setkey(env2,start,end)
  comparison=foverlaps(env1,env2,by.x=c('pos','pos_end'),by.y=c('start','end'),nomatch=NULL)


  hprobs=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',c,h))
  size=dim(hprobs[[1]])[1]
  n_markers=dim(hprobs[[1]])[2]

  comparison=comparison[comparison$marker %in% dimnames(hprobs[[1]])[[2]],]
  markers=comparison$marker
  # Get weights for null expectation
  # for each marker, count the number of founders in each haplotype


  all_weights=sapply(seq(1,n_markers), function(i) get_weighted_null(i,comparison))
  hsums=data.frame(matrix(unlist(lapply(hprobs, function(x) round(colSums(x)))),nrow=length(hprobs),byrow=T),stringsAsFactors=F)
  names(hsums)=markers


  p_chi=sapply(seq(1,n_markers), function(x) chisq.test(x=hsums[,x],p=all_weights[,x])$p.value)
#p_chi=t(p_chi)
  p_chi=data.frame(marker=markers,p_chi=unlist(p_chi),hagrp=h,stringsAsFactors=F)
  rownames(p_chi)=seq(1,dim(p_chi)[1])
  p_chi$pos=comparison$pos

  for(m in p_chi[-log10(p_chi$p_chi)>=bonf,]$marker){
    count=count+1
    index=which(colnames(hsums)==m)
    o=hsums[,index]
    e=all_weights[,index]*nrow(hprobs[[1]])
    a=data.frame(actual=o,expected=e,stringsAsFactors=F)
    haplo_counts[[count]]=list(marker=m,counts=a,haplotype_members=comparison[comparison$marker==m,..hap_founders],start=comparison[comparison$marker==m,]$start,end=comparison[comparison$marker==m,]$end)
  }

  all_p_chi=rbind(all_p_chi,p_chi)
  #total_sites=total_sites + (n_markers*h)
}

saveRDS(haplo_counts,sprintf('chr%s_observed_expected_founder_counts.rds',c))

names(all_p_chi)=c('marker','p_chi','hapgrp','pos')


png(sprintf('selection/haplotype_probs/bg%s_chi_squared_scan.png',c),width=1600, height=800)
print(ggplot(all_p_chi,aes(x=pos/1e6,y=-log10(p_chi))) + geom_point() + geom_hline(yintercept=bonf,color="red") + xlab("Position (Mb)") + ylab("Chi-Squared -log10(P-Value)") + ggtitle(sprintf("Chi-Squared Test for Representation of Haplotypes on Chromosome %s",c)))
dev.off()

# In genetic distance
all_p_chi$cM=gmap[match(all_p_chi$marker,gmap$marker)]$pos

png(sprintf('selection/haplotype_probs/bg%s_chi_squared_scan_cM.png',c),width=1600, height=800)
print(ggplot(all_p_chi,aes(x=cM,y=-log10(p_chi))) + geom_point() + geom_hline(yintercept=bonf,color="red") + xlab("Position (cM)") + ylab("Chi-Squared -log10(P-Value)") + ggtitle(sprintf("Chi-Squared Test for Representation of Haplotypes on Chromosome %s",c)))
dev.off()


fwrite(all_p_chi,sprintf('selection/haplotype_probs/bg%s_haplotype_chisq_results.txt',c),row.names=F,quote=F,sep='\t')
