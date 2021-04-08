#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
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

s_f_h_qtl<-function(name,label){
  sub=qtl[qtl$pheno_env_id == name,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method=="Founder_probs",]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=f_ses[[which(sapply(f_ses, function(x) x$qtl_id==name))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  rownames(fmelt)=founders
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')

  #Haplotype effect sizes
  line=sub[sub$Method=="Haplotype_probs",]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=hap_table[hap_table$SNP == line$highest_SNP, ]$HAPGRP
  hdata=h_ses[[which(sapply(h_ses, function(x) x$qtl_id==name))]]
  hsnp=hdata$snp
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==line$highest_SNP,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  hmelt$hapgrp=seq(1,h)
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  line=sub[sub$Method=="600K_SNP",]
  sdata=s_ses[[which(sapply(s_ses, function(x) x$qtl_id==name))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  smelt$allele=alleles
  names(smelt)=c('value','se','tvalue','variable','allele')

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  return(list(values=fmelt,id=name,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp,label=label))
}

s_only_qtl<-function(name,label){
  sub=qtl[qtl$pheno_env_id == name,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method=="600K_SNP",]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=s_notfh_ses[[which(sapply(s_notfh_ses, function(x) x$qtl_id==name & x$method=="Founder_probs"))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  rownames(fmelt)=founders
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')

  #Haplotype effect sizes
  #line=sub[sub$Method=="Haplotype_probs",]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  hdata=s_notfh_ses[[which(sapply(s_notfh_ses, function(x) x$qtl_id==name & x$method=="Haplotype_probs"))]]
  hsnp=hdata$snp
  h=hap_table[hap_table$SNP == hsnp, ]$HAPGRP
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==hsnp,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  rownames(hmelt)=seq(1,h)
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  #line=sub[sub$Method=="600K_SNP",]
  sdata=s_ses[[which(sapply(s_ses, function(x) x$qtl_id==name))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  smelt$allele=alleles
  names(smelt)=c('value','se','tvalue','variable','allele')

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  return(list(values=fmelt,id=name,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp,label=label))
}

f_only_qtl<-function(name,label){
  sub=qtl[qtl$pheno_env_id == name,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method=="Founder_probs",]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=f_ses[[which(sapply(f_ses, function(x) x$qtl_id==name))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  rownames(fmelt)=founders
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')

  #Haplotype effect sizes
  #line=sub[sub$Method=="Haplotype_probs",]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  hdata=f_notsh_ses[[which(sapply(f_notsh_ses, function(x) x$qtl_id==name & x$method=="Haplotype_probs"))]]
  hsnp=hdata$snp
  h=hap_table[hap_table$SNP == hsnp, ]$HAPGRP
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==hsnp,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  rownames(hmelt)=seq(1,h)
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  #line=sub[sub$Method=="600K_SNP",]
  sdata=f_notsh_ses[[which(sapply(f_notsh_ses, function(x) x$qtl_id==name & x$method=="600K_SNP"))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  #alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  #smelt$allele=alleles
  names(smelt)=c('value','se','tvalue','allele','variable')

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  return(list(values=fmelt,id=name,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp,label=label))
}

h_only_qtl<-function(name,label){
  sub=qtl[qtl$pheno_env_id == name,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method=="Haplotype_probs",]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=h_notsf_ses[[which(sapply(h_notsf_ses, function(x) x$qtl_id==name & x$method=="Founder_probs"))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  rownames(fmelt)=founders
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')

  #Haplotype effect sizes
  #line=sub[sub$Method=="Haplotype_probs",]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  hdata=h_ses[[which(sapply(h_ses, function(x) x$qtl_id==name))]]
  hsnp=hdata$snp
  h=hap_table[hap_table$SNP == hsnp, ]$HAPGRP
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==hsnp,]$pos
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  rownames(hmelt)=seq(1,h)
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  #line=sub[sub$Method=="600K_SNP",]
  sdata=h_notsf_ses[[which(sapply(h_notsf_ses, function(x) x$qtl_id==name & x$method=="600K_SNP"))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  #alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  #smelt$allele=alleles
  names(smelt)=c('value','se','tvalue','allele','variable')

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  return(list(values=fmelt,id=name,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp,label=label))
}

s_and_f_qtl<-function(name,label){
  sub=qtl[qtl$pheno_env_id == name,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method=="Founder_probs",]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=f_ses[[which(sapply(f_ses, function(x) x$qtl_id==name))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  rownames(fmelt)=founders
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')

  #Haplotype effect sizes
  line=sub[sub$Method=="Haplotype_probs",]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  hdata=noth_ses[[which(sapply(noth_ses, function(x) x$qtl_id==name))]]
  hsnp=hdata$snp
  h=hap_table[hap_table$SNP == hsnp, ]$HAPGRP
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==hsnp,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  rownames(hmelt)=seq(1,h)
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  line=sub[sub$Method=="600K_SNP",]
  sdata=s_ses[[which(sapply(s_ses, function(x) x$qtl_id==name))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  smelt$allele=alleles
  names(smelt)=c('value','se','tvalue','variable','allele')

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  return(list(values=fmelt,id=name,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp,label=label))
}

f_and_h_qtl<-function(name,label){
  sub=qtl[qtl$pheno_env_id == name,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method=="Founder_probs",]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=f_ses[[which(sapply(f_ses, function(x) x$qtl_id==name))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  rownames(fmelt)=founders
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')

  #Haplotype effect sizes
  line=sub[sub$Method=="Haplotype_probs",]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  hdata=h_ses[[which(sapply(h_ses, function(x) x$qtl_id==name))]]
  hsnp=hdata$snp
  h=hap_table[hap_table$SNP == hsnp, ]$HAPGRP
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==line$highest_SNP,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  rownames(hmelt)=seq(1,h)
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  line=sub[sub$Method=="600K_SNP",]
  sdata=nots_ses[[which(sapply(nots_ses, function(x) x$qtl_id==name))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  #alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  #smelt$allele=alleles
  names(smelt)=c('value','se','tvalue','allele','variable')


  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  return(list(values=fmelt,id=name,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp,label=label))
}


ids=unique(qtl_overlap$pheno_env_id)
count=1
plots=list()
for(i in seq(1,nrow(qtl_overlap))){
  name=qtl_overlap[i,]$pheno_env_id
  label=qtl_overlap[i,]$label
  if(label=='S_only'){
    plots[[count]]=s_only_qtl(name,label)
  }
  else if(label=='F_only'){
    plots[[count]]=f_only_qtl(name,label)
  }
  else if(label=="H_only"){
    plots[[count]]=h_only_qtl(name,label)
  }
  else if(label=="S_and_F"){
    plots[[count]]=s_and_f_qtl(name,label)
  }
  else if(label=="F_and_H"){
    plots[[count]]=f_and_h_qtl(name,label)
  }
  else{
    plots[[count]]=s_f_h_qtl(name,label)
  }
  count=count+1
}

saveRDS(plots,'GridLMM/effect_sizes/All_effect_sizes.rds')
