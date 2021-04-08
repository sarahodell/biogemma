#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
chr=as.numeric(args[[2]])
h=as.numeric(args[[3]])
cores=as.numeric(args[[4]])

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')

xfile=sprintf('../../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%.0f_filtered_haplogroup%.0f_probs.rds',chr,h)

if(file.exists(xfile)){
  X_list=readRDS(xfile)
  # Read in Kinship Matrix
  K=fread(sprintf('../K_matrices/K_matrix_chr%.0f.txt',chr),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)


  # Read in phenotypes
  # Grab the phenotype of interest and drop the genotypes not in the K matrix
  phenotypes=fread('../phenotypes_asi.csv',data.table=F)
  phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
  phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
  phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
  data=data[!is.na(data$y),]
  data$y=data$y - mean(data$y)

  m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
  data_blup = as.data.frame(ranef(m1)$ID2)
  data_blup$ID = rownames(data_blup)
  data_blup$y=data_blup$`(Intercept)`
  data_blup=data_blup[,c('ID','y')]
  inds=rownames(X_list[[1]])
  i=intersect(data_blup$ID,inds)
  # Read in the haplotype group probabilities
  # Filter genotypes that are not in the K matrix
  K=K[i,i]
  data_blup=data_blup[i,]
  # Run GridLMM
  null_model = GridLMM_ML(y~1+(1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F)

  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
  h2_start
  V_setup=null_model$setup

  Y=as.matrix(data_blup$y)
  X_cov=null_model$lmod$X
  dimx=dim(X_list[[1]])[2]
  X_list_ordered=lapply(X_list,function(x) array(x[i,],dim=c(length(i),dimx),dimnames=list(i,dimnames(X_list[[1]])[[2]])))

  X_list_null=NULL

  gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

  hinfo=data.frame(method="Haplotype_probs",phenotype=pheno,environment="ALL",chr=chr,h2=h2_start,hap=h,stringsAsFactors=F)
  fwrite(hinfo,'../heritabilities.txt',quote=F,sep='\t',row.names=F,append=T)

  saveRDS(gwas,sprintf('models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_ALL.rds',chr,h,pheno))
}else{
  print(sprintf("No markers at %s haplogroup %.0f",chr,h))
}

# Convert all very high and very low probabilities to 1 and 0, respectively
#X_list_full = lapply(X_list_ordered,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.99,1,ifelse(x[,i]<=1e-2,0,x[,i]))))
#for(i in 1:h){dimnames(X_list_full[[i]])[[2]]=dimnames(X_list_ordered[[i]])[[2]]}

#gwas_adjusted=gwas
#sums=lapply(X_list_full,function(x) colSums(x))
#for(i in 1:h){
#    s=sums[[i]]
#    t=dim(X_list_full[[i]])[1]-2
#    l=2
#    grab=which(s>t,s)
#    grab=c(grab,which(s<l,s))
#    grab=sort(grab)
#    beta=sprintf('beta.%.0f',seq(1,h))
#    gwas_adjusted[grab,beta]=0
#    print(grab)
#}

#saveRDS(gwas_adjusted,sprintf('models/Biogemma_chr%s_haplogrp%.0f_%s_x_ALL_adjusted.rds',chr,h,pheno))
