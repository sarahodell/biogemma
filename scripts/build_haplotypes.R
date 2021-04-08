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

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genofile=sprintf("genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds",c)
pmapfile=sprintf("genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv",c)
ibdfile=sprintf('ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',c)

# Read in IBD segments from get_ibd.R
ibd_segments=fread(ibdfile,data.table=T)
pmap=fread(pmapfile,data.table=F)
#founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


pr=readRDS(genofile)
dimnames(pr[[1]])[[2]]=founders
pmap=pmap[pmap$marker %in% dimnames(pr[[1]])[[3]],]
pr=lapply(pr, function(x) x[,,pmap$marker])
#pr=lapply(pr,function(x) x[,hap_founders,])

K=fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',c),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

phenotypes=fread('GridLMM/phenotypes_asi.csv',data.table=F)
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
pheno_IDs=unique(phenotypes$Genotype_code)

pr=lapply(pr,function(x) x[pheno_IDs,,])


size=dim(pr[[1]])[3]
samples=unlist(dimnames(pr[[1]])[[1]])

print(sprintf("Starting haplotype %s",h))
n_hap=h
n_ind=dim(pr[[1]])[1]
test=ibd_segments[ibd_segments$n_haps==n_hap,]
if(dim(test)[1]!=0){
  rownames(test)=seq(1,dim(test)[1])
  env1=pmap
  env1$pos_end=env1$pos
  env1=as.data.table(env1)
  env2=test
  env2$end=env2$end-1
  setkey(env2,start,end)
  comparison=foverlaps(env1,env2,by.x=c('pos','pos_end'),by.y=c('start','end'),nomatch=NULL)
  markers=comparison$marker
  hprobs=vector("list",length=n_hap)

  get_haploprobs<-function(i,comparison,pr){
    hprob=sapply(seq(1,n_ind), function(x) rowSums(t(pr[[1]][x,,comparison$marker]) * as.matrix(comparison[,..founders]==i)))
    colnames(hprob)=samples
    return(hprob)
  }
  for(i in 1:n_hap){
    hprobs[[i]]=t(get_haploprobs(i,comparison,pr))
  }


  outfile=sprintf("genotypes/probabilities/haplotype_probs/RefinedIBD_600K/raw/bg%s_refined_ibd_haplogroup%.0f_probs.rds",c,h)
  saveRDS(hprobs,outfile)
}else{
  print(sprintf("No markers in %s haplogroup %.0f",c,h))
}
