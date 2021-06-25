#!/usr/bin/env Rscript

library('data.table')
library('lme4')
library('lme4qtl')
library('emmeans')
library('ggplot2')
library('GridLMM')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="Founder_probs",]

rownames(qtl)=seq(1,nrow(qtl))
ses=list()
count=1
falsep=c()
freps=c()

for(i in 1:nrow(qtl)){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id
  cutoff=thresh[thresh$phenotyp==pheno & thresh$environment==env & thresh$method=="founder_probs",]$threshold
  #print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
  colnames(X) = founders
  rownames(X) = dimnames(founder_probs[[1]])[[1]]

  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    #rownames(phenotype)=seq(1,nrow(phenotype))
    #X = X[rownames(X) %in% phenotype$Genotype_code,]
    rownames(phenotype)=phenotype$Genotype_code
    i=intersect(phenotype$Genotype_code,rownames(X))
    X = X[i,]
    phenotype=phenotype[i,]
    subK=K[i,i]
    frep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    fkeep=founders[frep2>2]
    X=X[,fkeep]

    #GridLMM
    null_model = GridLMM_ML(y~1 + (1|Genotype_code),phenotype,relmat=list(Genotype_code=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m4))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'Genotype_code'

    Y=as.matrix(phenotype$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    if(-log10(gwas$p_value_ML)<cutoff & length(fkeep)<16){
      print(name)
      falsep=c(falsep,name)
      freps=rbind(freps,c(-log10(gwas$p_value_ML),cutoff,name,as.numeric(unlist(frep2))))
    }
  }else{
    phenotype = phenotype[!is.na(phenotype$y),]
    #phenotype$y=phenotype$y-mean(phenotype$y)
    m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m0)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    i=intersect(data_blup$ID,rownames(X))
    X = X[i,]
    data_blup=data_blup[i,]
    subK=K[i,i]
    frep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    fkeep=founders[frep2>2]
    X=X[,fkeep]

      #GRIDLMM
    null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m4))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'ID'


    Y=as.matrix(data_blup$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    if(-log10(gwas$p_value_ML)<cutoff & length(fkeep)<16){
      print(name)
      falsep=c(falsep,name)
      freps=rbind(freps,c(-log10(gwas$p_value_ML),cutoff,name,as.numeric(unlist(frep2))))
    }
  }
}

#Haplotypes
qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env=paste0(qtl$Phenotype,'_',qtl$Environment)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="Haplotype_probs",]
#qtl=qtl[qtl$Phenotype==pheno & qtl$Environment==env & qtl$Chromosome==as.numeric(chr),]
rownames(qtl)=seq(1,nrow(qtl))

falsep=c()
hreps=c()
for(i in 1:nrow(qtl)){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id
  #print(name)
  cutoff=thresh[thresh$phenotyp==pheno & thresh$environment==env & thresh$method=="haplotype_probs",]$threshold
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=results[results$SNP==snp,]$HAPGRP
  all_haps=paste0("HAPGRP_",seq(1,h))
  haplo_probs = readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(haplo_probs,function(x) x[,snp]))
  colnames(X) = paste0("HAPGRP_",seq(1,h))
  rownames(X) = dimnames(haplo_probs[[1]])[[1]]
  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=phenotype$Genotype_code
    i=intersect(phenotype$Genotype_code,rownames(X))
    X = X[i,]
    phenotype=phenotype[i,]
    subK=K[i,i]
    hrep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    hkeep=all_haps[hrep2>2]
    X=X[,hkeep]

    null_model = GridLMM_ML(y~1 + (1|Genotype_code),phenotype,relmat=list(Genotype_code=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m0))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'Genotype_code'

    Y=as.matrix(phenotype$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    if(-log10(gwas$p_value_ML)<cutoff & length(hkeep)<h){
      print(name)
      falsep=c(falsep,name)
      hreps=rbind(hreps,c(-log10(gwas$p_value_ML),cutoff,name,as.numeric(unlist(hrep2))))
    }
  }
  else{
    phenotype = phenotype[!is.na(phenotype$y),]
    #phenotype$y=phenotype$y-mean(phenotype$y)
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    i=intersect(data_blup$ID,rownames(X))
    X = X[i,]
    data_blup=data_blup[i,]
    subK=K[i,i]
    hrep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    hkeep=all_haps[hrep2>2]
    X=X[,hkeep]
    #GRIDLMM
    null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m0))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'ID'


    Y=as.matrix(data_blup$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    if(-log10(gwas$p_value_ML)<cutoff & length(hkeep)<h){
      print(name)
      falsep=c(falsep,name)
      hreps=rbind(hreps,c(-log10(gwas$p_value_ML),cutoff,name,as.numeric(unlist(hrep2))))
    }
  }
}
