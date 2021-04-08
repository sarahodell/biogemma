#!/usr/bin/env Rscript

library('data.table')
library('lme4')

K=fread('GridLMM/K_matrices/K_matrix_chr10.txt',data.table=F)
inds=K[,1]
qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
#qtl=qtl$


phenotypes=fread('GridLMM/phenotypes_asi.csv',data.table=F)
#phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% inds,]
#data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
ft_days=c('male_flowering_days','female_flowering_days')
qtl_overlap=qtl_overlap[!(qtl_overlap$Phenotype %in% ft_days),]
rownames(qtl_overlap)=seq(1,nrow(qtl_overlap))
#qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)

f_ses=readRDS('GridLMM/effect_sizes/Founder_prob_QTL_SEs.rds')
h_ses=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
s_ses=readRDS('GridLMM/effect_sizes/600K_SNP_QTL_SEs.rds')

#S only
s_notfh_ses=readRDS('GridLMM/effect_sizes/S_only_not_F_and_H_QTL_SEs.rds')
s_notfh=fread('GridLMM/result_tables/SNP_only_not_F_and_H_highest_peaks.txt',data.table=F)

# F only
f_notsh_ses=readRDS('GridLMM/effect_sizes/F_only_not_S_and_H_QTL_SEs.rds')
f_notsh=fread('GridLMM/result_tables/F_only_not_S_and_H_highest_peaks.txt',data.table=F)

# H only
h_notsf_ses=readRDS('GridLMM/effect_sizes/H_only_not_S_and_F_QTL_SEs.rds')
h_notsf=fread('GridLMM/result_tables/H_only_not_S_and_F_highest_peaks.txt',data.table=F)

# S and F
noth_ses=readRDS('GridLMM/effect_sizes/S_and_F_not_H_QTL_SEs.rds')
noth=fread('GridLMM/result_tables/S_and_F_not_H_highest_peaks.txt')

nots_ses=readRDS('GridLMM/effect_sizes/F_and_H_not_S_QTL_SEs.rds')
nots=fread('GridLMM/result_tables/F_and_H_not_S_highest_peaks.txt')

s_f_h_ids=qtl_overlap[qtl_overlap$label=="S_F_H",]$pheno_env_id
s_only_ids=qtl_overlap[qtl_overlap$label=="S_only",]$pheno_env_id
f_only_ids=qtl_overlap[qtl_overlap$label=="F_only",]$pheno_env_id
h_only_ids=qtl_overlap[qtl_overlap$label=="H_only",]$pheno_env_id
s_and_f_ids=qtl_overlap[qtl_overlap$label=="S_and_F",]$pheno_env_id
f_and_h_ids=qtl_overlap[qtl_overlap$label=="F_and_H",]$pheno_env_id


#count=1
sqtl=qtl[qtl$Method=="600K_SNP",]
rownames(sqtl)=seq(1,nrow(sqtl))
for(i in seq(1,nrow(sqtl))){
  line=sqtl[i,]
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  pe_id=line$pheno_env_id
  highest_SNP=line$highest_SNP
  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
  if(env=="ALL"){
    m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
    data_blup = as.data.frame(ranef(m1)$ID2)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    tot_var=var(data_blup$y)
  }else{
    data=data[data$Loc.Year.Treat==env,]
    data=data[!is.na(data$y),]
    tot_var=var(data$y)
  }
  #smodel=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%s_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
  #smodel=smodel$results
  sgeno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
  #ses=readRDS('GridLMM/effect_sizes/600K_SNP_QTL_SEs.rds')
  ses=s_ses[[which(unlist(lapply(s_ses,function(x) x$qtl_id==pe_id)))]]
  sgeno=sgeno[,c('ind',highest_SNP)]
  sgeno=sgeno[sgeno$ind %in% inds,]
  #N=nrow(sgeno)
  #MAF=sum(sgeno[,2])/N
  #if(MAF > 0.5){
  #  MAF=1-MAF
  #}
  #smodel=smodel[smodel$X_ID==highest_SNP,]
  betas=ses$SE
  beta=betas[2,]$value
  se_beta=betas[2,]$se
  snp_var=var(sgeno[,highest_SNP])
  pve=(beta^2 * snp_var)/tot_var
  #pve=(2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF))
  pves=rbind(pves,c('600K_SNP',pe_id,(beta^2 * snp_var),tot_var,pve,T))
}

fqtl=qtl[qtl$Method=="Founder_probs",]
rownames(fqtl)=seq(1,nrow(fqtl))
for(i in seq(1,nrow(fqtl))){
  line=fqtl[i,]
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  pe_id=line$pheno_env_id
  highest_SNP=line$highest_SNP
  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

  if(env=="ALL"){
    m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
    data_blup = as.data.frame(ranef(m1)$ID2)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    tot_var=var(data_blup$y)
  }else{
    data=data[data$Loc.Year.Treat==env,]
    data=data[!is.na(data$y),]
    tot_var=var(data$y)
  }
  #smodel=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%s_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
  #smodel=smodel$results
  fgeno=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  #fes=readRDS('GridLMM/effect_sizes/Founder_prob_QTL_SEs.rds')
  fes=f_ses[[which(unlist(lapply(f_ses,function(x) x$qtl_id==pe_id)))]]
  fgeno=data.frame(lapply(fgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
  fgeno=fgeno[rownames(fgeno) %in% inds,]
  N=nrow(fgeno)
  #MAF=sum(fgeno[,2])/N
  #if(MAF > 0.5){
  #  MAF=1-MAF
  #}
  #smodel=smodel[smodel$X_ID==highest_SNP,]
  betas=fes$SE
  rownames(betas)=names(fgeno)
  betas=betas[complete.cases(betas),]
  founders=rownames(betas)
  est=sapply(seq(1,N),function(x) sum(fgeno[x,founders] * betas$value))
  pve=var(est)/tot_var
  #se_beta=betas[2,]$se
  #snp_var=var(sgeno[,highest_SNP])
  #pve=(2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF))
  pves=rbind(pves,c('Founder_probs',pe_id,var(est),tot_var,pve,T))
}

hqtl=qtl[qtl$Method=="Haplotype_probs",]
rownames(hqtl)=seq(1,nrow(hqtl))
for(i in seq(1,nrow(hqtl))){
  line=hqtl[i,]
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  pe_id=line$pheno_env_id
  highest_SNP=line$highest_SNP
  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
  if(env=="ALL"){
    m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
    data_blup = as.data.frame(ranef(m1)$ID2)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    tot_var=var(data_blup$y)
  }else{
    data=data[data$Loc.Year.Treat==env,]
    data=data[!is.na(data$y),]
    tot_var=var(data$y)
  }
  #smodel=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%s_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
  #smodel=smodel$results
  hresults=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env))
  h=hresults[hresults$SNP==highest_SNP,]$HAPGRP
  hgeno=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
  hes=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
  hes=hes[[which(unlist(lapply(hes,function(x) x$qtl_id==pe_id)))]]
  hgeno=data.frame(lapply(hgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
  names(hgeno)=paste0('HAPGRP_',seq(1,h))
  hgeno=hgeno[rownames(hgeno) %in% inds,]
  N=nrow(hgeno)
  #MAF=sum(fgeno[,2])/N
  #if(MAF > 0.5){
  #  MAF=1-MAF
  #}
  #smodel=smodel[smodel$X_ID==highest_SNP,]
  betas=hes$SE
  #rownames(betas)=names(hgeno)
  betas=betas[complete.cases(betas),]
  hapgrps=rownames(betas)
  est=sapply(seq(1,N),function(x) sum(hgeno[x,hapgrps] * betas$value))
  pve=var(est)/tot_var
  #se_beta=betas[2,]$se
  #snp_var=var(sgeno[,highest_SNP])
  #pve=(2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF))
  pves=rbind(pves,c('Haplotype_probs',pe_id,var(est),tot_var,pve,T))
}
pves=as.data.frame(pves,stringsAsFactors=F)
names(pves)=c('method','pheno_env_id','site_var','tot_pheno_var','pve','sig')



# Non-QTL variance explained

# S only, not F and H
for(i in seq(1,length(s_notfh_ses))){
  print(i)
  fes=s_notfh_ses[[i]]
  line=s_notfh[i,]
  method=line$method
  if(method=="Founder_probs"){
    pheno=line$phenotype
    env=line$environment
    chr=as.character(line$chrom)
    pe_id=line$pheno_env_id
    highest_SNP=line$highest_SNP
    data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
    if(env=="ALL"){
      m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
      data_blup = as.data.frame(ranef(m1)$ID2)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      tot_var=var(data_blup$y)
    }
    else{
      data=data[data$Loc.Year.Treat==env,]
      data=data[!is.na(data$y),]
      tot_var=var(data$y)
    }
    fgeno=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
    fgeno=data.frame(lapply(fgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
    fgeno=fgeno[rownames(fgeno) %in% inds,]
    N=nrow(fgeno)
    betas=fes$SE
    rownames(betas)=names(fgeno)
    betas=betas[complete.cases(betas),]
    founders=rownames(betas)
    est=sapply(seq(1,N),function(x) sum(fgeno[x,founders] * betas$value))
    pve=var(est)/tot_var
    #pve=(2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF))
    pves=rbind(pves,c('Founder_probs',pe_id,var(est),tot_var,pve,F))
  }
  else{
    pheno=line$phenotype
    env=line$environment
    chr=as.character(line$chrom)
    pe_id=line$pheno_env_id
    highest_SNP=line$highest_SNP
    data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
    if(env=="ALL"){
      m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
      data_blup = as.data.frame(ranef(m1)$ID2)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      tot_var=var(data_blup$y)
    }
    else{
      data=data[data$Loc.Year.Treat==env,]
      data=data[!is.na(data$y),]
      tot_var=var(data$y)
    }
    hresults=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env))
    h=hresults[hresults$SNP==highest_SNP,]$HAPGRP
    hgeno=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
    hgeno=data.frame(lapply(hgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
    names(hgeno)=paste0('HAPGRP_',seq(1,h))
    hgeno=hgeno[rownames(hgeno) %in% inds,]
    N=nrow(hgeno)
    betas=fes$SE
    #rownames(betas)=names(hgeno)
    betas=betas[complete.cases(betas),]
    hapgrps=rownames(betas)
    est=sapply(seq(1,N),function(x) sum(hgeno[x,hapgrps] * betas$value))
    pve=var(est)/tot_var
    #pve=(2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF))
    pves=rbind(pves,c('Haplotype_probs',pe_id,var(est),tot_var,pve,F))
  }
}

# F only, not S and H
for(i in seq(1,length(f_notsh_ses))){
  fes=f_notsh_ses[[i]]
  line=f_notsh[i,]
  method=line$method
  if(method=="600K_SNP"){
    pheno=line$phenotype
    env=line$environment
    chr=as.character(line$chrom)
    pe_id=line$pheno_env_id
    highest_SNP=line$highest_SNP
    data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
    if(env=="ALL"){
      m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
      data_blup = as.data.frame(ranef(m1)$ID2)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      tot_var=var(data_blup$y)
    }else{
      data=data[data$Loc.Year.Treat==env,]
      data=data[!is.na(data$y),]
      tot_var=var(data$y)
    }
    sgeno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
    #ses=s_ses[[which(unlist(lapply(s_ses,function(x) x$qtl_id==pe_id)))]]
    sgeno=sgeno[,c('ind',highest_SNP)]
    sgeno=sgeno[sgeno$ind %in% inds,]
    betas=fes$SE
    beta=betas[2,]$value
    se_beta=betas[2,]$se
    snp_var=var(sgeno[,highest_SNP])
    pve=(beta^2 * snp_var)/tot_var
    pves=rbind(pves,c('600K_SNP',pe_id,(beta^2 * snp_var),tot_var,pve,F))
  }
  else{
    pheno=line$phenotype
    env=line$environment
    chr=as.character(line$chrom)
    pe_id=line$pheno_env_id
    highest_SNP=line$highest_SNP
    data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

    if(env=="ALL"){
      m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
      data_blup = as.data.frame(ranef(m1)$ID2)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      tot_var=var(data_blup$y)
    }else{
      data=data[data$Loc.Year.Treat==env,]
      data=data[!is.na(data$y),]
      tot_var=var(data$y)
    }
    hresults=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env))
    h=hresults[hresults$SNP==highest_SNP,]$HAPGRP
    hgeno=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
    hgeno=data.frame(lapply(hgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
    names(hgeno)=paste0('HAPGRP_',seq(1,h))
    hgeno=hgeno[rownames(hgeno) %in% inds,]
    N=nrow(hgeno)
    betas=fes$SE
    betas=betas[complete.cases(betas),]
    hapgrps=rownames(betas)
    est=sapply(seq(1,N),function(x) sum(hgeno[x,hapgrps] * betas$value))
    pve=var(est)/tot_var
    pves=rbind(pves,c('Haplotype_probs',pe_id,var(est),tot_var,pve,F))
  }
}

# H only, not F and S
for(i in seq(1,length(h_notsf_ses))){
  fes=h_notsf_ses[[i]]
  line=h_notsf[i,]
  method=line$method
  if(method=="600K_SNP"){
    pheno=line$phenotype
    env=line$environment
    chr=as.character(line$chrom)
    pe_id=line$pheno_env_id
    highest_SNP=line$highest_SNP
    data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

    if(env=="ALL"){
      m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
      data_blup = as.data.frame(ranef(m1)$ID2)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      tot_var=var(data_blup$y)
    }else{
      data=data[data$Loc.Year.Treat==env,]
      data=data[!is.na(data$y),]
      tot_var=var(data$y)
    }
    sgeno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
    #ses=s_ses[[which(unlist(lapply(s_ses,function(x) x$qtl_id==pe_id)))]]
    sgeno=sgeno[,c('ind',highest_SNP)]
    sgeno=sgeno[sgeno$ind %in% inds,]
    betas=fes$SE
    beta=betas[2,]$value
    se_beta=betas[2,]$se
    snp_var=var(sgeno[,highest_SNP])
    pve=(beta^2 * snp_var)/tot_var
    pves=rbind(pves,c('600K_SNP',pe_id,(beta^2 * snp_var),tot_var,pve,F))
  }
  else{
    pheno=line$phenotype
    env=line$environment
    chr=as.character(line$chrom)
    pe_id=line$pheno_env_id
    highest_SNP=line$highest_SNP
    data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

    if(env=="ALL"){
      m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
      data_blup = as.data.frame(ranef(m1)$ID2)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      tot_var=var(data_blup$y)
    }else{
      data=data[data$Loc.Year.Treat==env,]
      data=data[!is.na(data$y),]
      tot_var=var(data$y)
    }
    fgeno=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
    fgeno=data.frame(lapply(fgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
    fgeno=fgeno[rownames(fgeno) %in% inds,]
    N=nrow(fgeno)
    betas=fes$SE
    rownames(betas)=names(fgeno)
    betas=betas[complete.cases(betas),]
    founders=rownames(betas)
    est=sapply(seq(1,N),function(x) sum(fgeno[x,founders] * betas$value))
    pve=var(est)/tot_var
    #pve=(2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF))
    pves=rbind(pves,c('Founder_probs',pe_id,var(est),tot_var,pve,F))
  }
}

# S and F, not H
for(i in seq(1,length(noth_ses))){
  fes=noth_ses[[i]]
  line=noth[i,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  pe_id=line$pheno_env_id
  highest_SNP=line$highest_SNP
  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
  if(env=="ALL"){
    m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
    data_blup = as.data.frame(ranef(m1)$ID2)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    tot_var=var(data_blup$y)
  }else{
    data=data[data$Loc.Year.Treat==env,]
    data=data[!is.na(data$y),]
    tot_var=var(data$y)
  }
  hresults=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env))
  h=hresults[hresults$SNP==highest_SNP,]$HAPGRP
  hgeno=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
  #hes=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
  #hes=hes[[which(unlist(lapply(hes,function(x) x$qtl_id==pe_id)))]]
  hgeno=data.frame(lapply(hgeno,function(x) x[,highest_SNP]),stringsAsFactors=F)
  names(hgeno)=paste0('HAPGRP_',seq(1,h))
  hgeno=hgeno[rownames(hgeno) %in% inds,]
  N=nrow(hgeno)
  betas=fes$SE
  #rownames(betas)=names(hgeno)
  betas=betas[complete.cases(betas),]
  hapgrps=rownames(betas)
  est=sapply(seq(1,N),function(x) sum(hgeno[x,hapgrps] * betas$value))
  pve=var(est)/tot_var
  pves=rbind(pves,c('Haplotype_probs',pe_id,var(est),tot_var,pve,F))
}

# F and H, not S
for(i in seq(1,length(nots_ses))){
  fes=nots_ses[[i]]
  line=nots[i,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  pe_id=line$pheno_env_id
  highest_SNP=line$highest_SNP
  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

  if(env=="ALL"){
    m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
    data_blup = as.data.frame(ranef(m1)$ID2)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    tot_var=var(data_blup$y)
  }else{
    data=data[data$Loc.Year.Treat==env,]
    data=data[!is.na(data$y),]
    tot_var=var(data$y)
  }
  sgeno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
  #ses=s_ses[[which(unlist(lapply(s_ses,function(x) x$qtl_id==pe_id)))]]
  sgeno=sgeno[,c('ind',highest_SNP)]
  sgeno=sgeno[sgeno$ind %in% inds,]
  betas=fes$SE
  beta=betas[2,]$value
  se_beta=betas[2,]$se
  snp_var=var(sgeno[,highest_SNP])
  pve=(beta^2 * snp_var)/tot_var
  pves=rbind(pves,c('600K_SNP',pe_id,(beta^2 * snp_var),tot_var,pve,F))
}

fwrite(pves,'GridLMM/Biogemma_QTL_variance_explained.txt',row.names=F,quote=F,sep='\t')
