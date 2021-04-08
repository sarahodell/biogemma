#!/usr/bin/env Rscript

library('data.table')
library('lme4')
library('lme4qtl')
library('emmeans')
library('ggplot2')
#library('GridLMM')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)

qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
qtl$grouptype=qtl_overlap[match(qtl$pheno_env_id,qtl_overlap$pheno_env_id),]$label

s_notfh=fread('GridLMM/result_tables/SNP_only_not_F_and_H_highest_peaks.txt',data.table=F)
f_notsh=fread('GridLMM/result_tables/F_only_not_S_and_H_highest_peaks.txt',data.table=F)
h_notsf=fread('GridLMM/result_tables/H_only_not_S_and_F_highest_peaks.txt',data.table=F)
noth=fread('GridLMM/result_tables/S_and_F_not_H_highest_peaks.txt')
nots=fread('GridLMM/result_tables/F_and_H_not_S_highest_peaks.txt')

snp_se=function(i,focal,qtl){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  X = fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
  rownames(X)=X$ind
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  test=X[,names(X)[2],drop=F]
  X = X[,snp,drop=F]

  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=phenotype$Genotype_code
    i=intersect(phenotype$Genotype_code,rownames(X))
    X = X[i,,drop=F]
    phenotype=phenotype[i,]
    subK=K[i,i]
    m0 = relmatLmer(y ~ 1 + unlist(X) + (1|Genotype_code),data=phenotype,
    relmat = list(Genotype_code=subK))
    se=as.data.frame(summary(m0)$coef)
    names(se)=c('value','se','tvalue')
    se[2,]$value=se[-1,]$value + se[1,]$value
    se$variable_f=factor(c(0,1),levels=c(0,1))
    alt1=se$value[-1]+se$value[1]
  }
  else{
    phenotype = phenotype[!is.na(phenotype$y),]
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]

    i=intersect(data_blup$ID,rownames(X))
    X = X[i,,drop=F]
    data_blup=data_blup[i,]
    subK=K[i,i]
    m0 = relmatLmer(y ~ 1 + unlist(X) + (1|ID),data=data_blup,relmat = list(ID=subK))
    se=as.data.frame(summary(m0)$coef)
    names(se)=c('value','se','tvalue')
    se[2,]$value=se[-1,]$value + se[1,]$value
    se$variable_f=factor(c(0,1),levels=c(0,1))
    alt1=se$value[-1]+se$value[1]
  }
  return(list(SE=se,qtl_id=name,snp=snp,keptsnp=snp,pos=pos,method="600K_SNP",focal=focal))
}

founder_se=function(i,focal,qtl){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  lowrep=readRDS(sprintf('genotypes/probabilities/geno_probs/dropped/bg%s_dropped_markers_genoprobs.rds',chr))
  if(snp %in% dimnames(founder_probs[[1]])[[2]]){
    keptm=snp
  }else{
    findm=which(unlist(lapply(lowrep, function(x) snp %in% x$linked)))
    if(length(findm)==0){
      keptm=snp
    }else{
      keptm=lowrep[[findm]]$marker
    }
  }
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(founder_probs,function(x) x[,keptm]))
  colnames(X) = founders
  rownames(X) = dimnames(founder_probs[[1]])[[1]]
  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=phenotype$Genotype_code
    i=intersect(phenotype$Genotype_code,rownames(X))
    X = X[i,]
    phenotype=phenotype[i,]
    subK=K[i,i]
    frep=colSums(X)
    fkeep=founders[frep>1]
    X=X[,fkeep]
    m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
    se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
    names(se4)=c('value','se','tvalue')
    rownames(se4)=fkeep
    se4$founder=rownames(se4)
    se4$variable_f=factor(se4$founder,levels=se4$founder)
  }
  else{
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
    frep=colSums(X)
    fkeep=founders[frep>1]
    X=X[,fkeep]
    m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
    se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
    rownames(se4)=fkeep
    names(se4)=c('value','se','tvalue')
    se4$founder=rownames(se4)
    se4$variable_f=factor(se4$founder,levels=se4$founder)
  }
  return(list(SE=se4,qtl_id=name,snp=snp,keptsnp=keptm,pos=pos,method="Founder_prob",focal=focal))
}

haplotype_se=function(i,focal,qtl){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  base_list=c(7,7,8,6,7,6,7,7,6,8)
  base=base_list[as.numeric(chr)]
  for(h in base:16){
    hapfile=sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h)
    if(file.exists(hapfile)){
      haplo_probs = readRDS(hapfile)
      if(snp %in% dimnames(haplo_probs[[1]])[[2]]){
        hapgrp=h
        keptm=snp
      }else{
        lowrep=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/dropped/bg%s_haplogroup%.0f_dropped_markers.rds',chr,h))
        findm=which(unlist(lapply(lowrep, function(x) snp %in% x$linked)))
        if(length(findm)!=0){
          print(findm)
          keptm=lowrep[[findm]]$marker
          hapgrp=h
        }
      }
    }
  }
  hapfile=sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,hapgrp)
  haplo_probs = readRDS(hapfile)
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(haplo_probs,function(x) x[,keptm]))
  allhaps= paste0("HAPGRP_",seq(1,hapgrp))
  colnames(X) = allhaps
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
    hreps=colSums(X)
    hkeep=allhaps[hreps>1]
    X=X[,hkeep]
    m0 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
    se4=as.data.frame(summary(m0)$coef,stringsAsFactors=F)
    names(se4)=c('value','se','tvalue')
    rownames(se4)=hkeep
    se4$hapgrp=rownames(se4)
    se4$variable_f=factor(se4$hapgrp,levels=hkeep)
  }
  else{
    phenotype = phenotype[!is.na(phenotype$y),]
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    i=intersect(data_blup$ID,rownames(X))
    X = X[i,]
    data_blup=data_blup[i,]
    subK=K[i,i]
    hreps=colSums(X)
    hkeep=allhaps[hreps>1]
    X=X[,hkeep]
    m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK))
    se4=as.data.frame(summary(m0)$coef,stringsAsFactors=F)
    rownames(se4)=hkeep
    names(se4)=c('value','se','tvalue')
    se4$hapgrp=rownames(se4)
    se4$variable_f=factor(se4$hapgrp,levels=hkeep)
  }
  list(SE=se4,qtl_id=name,snp=snp,keptsnp=keptm,pos=pos,method="Haplotype_probs",hapgrp=hapgrp,focal=focal)
}

ses=list()
count=1
for(i in 1:nrow(qtl)){
  line=qtl[i,]
  method=line$Method
  if(method=="600K_SNP"){
    res=snp_se(i,"600K_SNP",qtl)
    res2=founder_se(i,"600K_SNP",qtl)
    res3=haplotype_se(i,"600K_SNP",qtl)
  }else if(method=="Founder_probs"){
    res=snp_se(i,"Founder_probs",qtl)
    res2=founder_se(i,"Founder_probs",qtl)
    res3=haplotype_se(i,"Founder_probs",qtl)
  }else{
    res=snp_se(i,"Haplotype_probs",qtl)
    res2=founder_se(i,"Haplotype_probs",qtl)
    res3=haplotype_se(i,"Haplotype_probs",qtl)
  }
  ses[[count]]=res
  count=count+1
  ses[[count]]=res2
  count=count+1
  ses[[count]]=res3
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/All_QTL_ES.rds')
