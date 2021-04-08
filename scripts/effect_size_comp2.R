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
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)

f_ses=readRDS('GridLMM/effect_sizes/Founder_prob_QTL_SEs.rds')
h_ses=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
s_ses=readRDS('GridLMM/effect_sizes/600K_SNP_QTL_SEs.rds')

#S only
s_notfh_ses=readRDS('GridLMM/effect_sizes/founder_ES/S_only_not_F_and_H_QTL_SEs.rds')
s_notfh=fread('GridLMM/result_tables/SNP_only_not_F_and_H_highest_peaks.txt',data.table=F)

# F only
f_notsh_ses=readRDS('GridLMM/effect_sizes/founder_ES/F_only_not_S_and_H_QTL_SEs.rds')
f_notsh=fread('GridLMM/result_tables/F_only_not_S_and_H_highest_peaks.txt',data.table=F)

# H only
h_notsf_ses=readRDS('GridLMM/effect_sizes/founder_ES/H_only_not_S_and_F_QTL_SEs.rds')
h_notsf=fread('GridLMM/result_tables/H_only_not_S_and_F_highest_peaks.txt',data.table=F)

# S and F
noth_ses=readRDS('GridLMM/effect_sizes/founder_ES/S_and_F_not_H_QTL_SEs.rds')
noth=fread('GridLMM/result_tables/S_and_F_not_H_highest_peaks.txt')

nots_sesreadRDS('GridLMM/effect_sizes/founder_ES/F_and_H_not_S_QTL_SEs.rds')
nots=fread('GridLMM/result_tables/F_and_H_not_S_highest_peaks.txt')

s_f_h_ids=qtl_overlap[qtl_overlap$label=="S_F_H",]$pheno_env_id
s_only_ids=qtl_overlap[qtl_overlap$label=="S_only",]$pheno_env_id
f_only_ids=qtl_overlap[qtl_overlap$label=="F_only",]$pheno_env_id
h_only_ids=qtl_overlap[qtl_overlap$label=="H_only",]$pheno_env_id
s_and_f_ids=qtl_overlap[qtl_overlap$label=="S_and_F",]$pheno_env_id
f_and_h_ids=qtl_overlap[qtl_overlap$label=="F_and_H",]$pheno_env_id



nots_plots=list()
count=0

# In founders but not SNPs
m1=c('Founder_probs','600K_SNP')
nots_ids=unique(f_nots$pheno_env_id)
for(q in f_nots_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method==m1[1],]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  #fresults=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
  #feffects=fresults[fresults$X_ID==line$highest_SNP,c(2,6:21)]
  #names(feffects)=c('X_ID',founders)
  #fmelt=melt(feffects,'X_ID')


  snp=f_nots[f_nots$pheno_env_id==q,]$highest_SNP
  sresults=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%.0f_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
  seffects=sresults$results[sresults$results$X_ID==line$highest_SNP,c(2,6:7)]
  names(seffects)=c('X_ID',seq(0,1))
  smelt=melt(seffects,'X_ID')
  smelt$method="600K_SNP"
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f_121718.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[founders,line$highest_SNP]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$founders=founders
  smelt$variable=as.character(smelt$variable)
  names(smelt)=c('X_ID','allele','value','method','variable')
  smelt=smelt[,c('X_ID','variable','value','method','allele')]
  smelt=smelt[order(smelt$allele),]
  rownames(smelt)=seq(1,nrow(smelt))
  fmelt=as.data.frame(f_ses[[which(sapply(f_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
  fmelt$variable=founders
  fmelt$method="Founder_probs"
  #fmelt$allele=as.character(seq(3:18))
  names(fmelt)=c('value','se','tvalue','variable','method')
  fmelt$allele=smelt$allele
  fmelt$allele_f=factor(paste0(fmelt$variable,'_',fmelt$allele),levels=paste0(fmelt$variable,'_',fmelt$allele))
  #allmelt=rbind(smelt,fmelt)
  #smelt$f_value=fmelt$value
  count=count+1

  f_nots_plots[[count]]=ggplot(fmelt,aes(x=allele_f,y=value)) + geom_bar(aes(fill=allele),stat="identity") +  geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+ ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Founder effect Size of %s (QTL in Founder not SNP)',q)) + theme(axis.text.x=element_text(size=5),title=element_text(size=10))
  #f_nots_plots[[count]]=ggplot(allmelt,aes(x=variable,y=value)) + geom_bar(aes(fill=allele),stat="identity")+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Highest SNP effect Size of %s (QTL in Founder not SNP)',q))

}

pdf('GridLMM/effect_sizes/Founder_notSNP_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(f_nots_plots)){
  print(f_nots_plots[[i]])
}
dev.off()

# In SNPS not in Founders
s_notf_ids=unique(s_notf$pheno_env_id)
s_notf_plots=list()
count=0
for(q in s_notf_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method==m1[2],]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  snp=s_notf[s_notf$pheno_env_id==q,]$highest_SNP
  fmelt=as.data.frame(s_notf_ses[[which(sapply(s_notf_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)

  #fresults=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
  #feffects=fresults[fresults$X_ID==snp,c(2,6:21)]
  #names(feffects)=c('X_ID',founders)
  #fmelt=melt(feffects,'X_ID')
  fmelt$variable=founders
  fmelt$method="Founder_probs"
  fmelt$allele=as.character(seq(3:18))
  names(fmelt)=c('value','se','tvalue','variable','method')

  sresults=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%.0f_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
  seffects=sresults$results[sresults$results$X_ID==line$highest_SNP,c(2,6:7)]
  names(seffects)=c('X_ID',seq(0,1))
  smelt=melt(seffects,'X_ID')
  smelt$method="600K_SNP"
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f_121718.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[founders,line$highest_SNP]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$founders=founders
  smelt$variable=as.character(smelt$variable)
  names(smelt)=c('X_ID','allele','value','method','variable')
  smelt=smelt[,c('X_ID','variable','value','method','allele')]
  smelt=smelt[order(smelt$allele),]
  rownames(smelt)=seq(1,nrow(smelt))
  fmelt$allele=smelt$allele
  fmelt$allele_f=factor(paste0(fmelt$variable,'_',fmelt$allele),levels=paste0(fmelt$variable,'_',fmelt$allele))
  #allmelt=rbind(smelt,fmelt)
  #smelt$f_value=fmelt$value
  count=count+1

  s_notf_plots[[count]]=ggplot(fmelt,aes(x=allele_f,y=value)) + geom_bar(aes(fill=allele),stat="identity") + geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+ ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Founder effect Size of %s (QTL in SNP not Founder)',q)) + theme(axis.text.x=element_text(size=5),title=element_text(size=10))
  #f_nots_plots[[count]]=ggplot(allmelt,aes(x=variable,y=value)) + geom_bar(aes(fill=allele),stat="identity")+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Highest SNP effect Size of %s (QTL in Founder not SNP)',q))

}

pdf('GridLMM/effect_sizes/SNP_notFounders_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(s_notf_plots)){
  print(s_notf_plots[[i]])
}
dev.off()

f_noth_ids=unique(f_noth$pheno_env_id)
f_noth_plots=list()
count=0
m2=c('Founder_probs','Haplotype_probs')
for(q in f_noth_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method==m2[1],]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fmelt=as.data.frame(f_ses[[which(sapply(f_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
  fmelt$variable=founders
  rownames(fmelt)=founders
  #fresults=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
  #feffects=fresults[fresults$X_ID==line$highest_SNP,c(2,6:21)]
  #names(feffects)=c('X_ID',founders)
  #fmelt=melt(feffects,'X_ID')
  fmelt$method="Founder_probs"
  fmelt$hapgrp=factor(seq(1,16),order(seq(1,16)))

  snp=f_noth[f_noth$pheno_env_id==q,]$highest_SNP
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=hap_table[hap_table$SNP == snp, ]$HAPGRP
  hresults=readRDS(sprintf('GridLMM/GridLMM_haplotypes/models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',chr,h,pheno,env))
  heffects=hresults[hresults$X_ID==snp,c(2,6:(h+5))]
  names(heffects)=c("X_ID",seq(1,h))
  hmelt=melt(heffects,"X_ID")
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$founder=founders
  names(hmelt)=c('X_ID','hapgrp','value','variable')
  hmelt$hapgrp=factor(hmelt$hapgrp,order(seq(1,h)))
  hmelt$method="Haplotype_probs"
  hmelt=hmelt[,c('X_ID','variable','value','method','hapgrp')]
  fmelt$hapgrp=factor(hmelt$hapgrp,levels=seq(1,h))
  names(fmelt)=c('value','se','tvalue','variable','method','hapgrp')
  fmelt$variable_f=factor(fmelt$variable,levels=fmelt[order(fmelt$hapgrp),]$variable)
  count=count+1

  f_noth_plots[[count]]=ggplot(fmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity") + geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge()) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Founder effect Size of %s (QTL in Founder not Haplotype)',q)) + theme(axis.text.x=element_text(size=5),title=element_text(size=10))
  #f_nots_plots[[count]]=ggplot(allmelt,aes(x=variable,y=value)) + geom_bar(aes(fill=allele),stat="identity")+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Highest SNP effect Size of %s (QTL in Founder not SNP)',q))

}

pdf('GridLMM/effect_sizes/Founder_notHaplotype_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(f_noth_plots)){
  print(f_noth_plots[[i]])
}
dev.off()


h_notf_ids=unique(h_notf$pheno_env_id)
h_notf_plots=list()
count=0
m2=c('Founder_probs','Haplotype_probs')
for(q in h_notf_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method==m2[2],]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  snp=h_notf[h_notf$pheno_env_id==q,]$highest_SNP
  fmelt=as.data.frame(h_notf_ses[[which(sapply(h_notf_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
  fmelt$variable=founders
  #fresults=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
  #feffects=fresults[fresults$X_ID==snp,c(2,6:21)]
  #names(feffects)=c('X_ID',founders)
  #fmelt=melt(feffects,'X_ID')
  fmelt$method="Founder_probs"
  #fmelt$hapgrp=factor(seq(1,16),order(seq(1,16)))
  hmelt=as.data.frame(h_ses[[which(sapply(h_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)


  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=hap_table[hap_table$SNP == line$highest_SNP, ]$HAPGRP
  #hresults=readRDS(sprintf('GridLMM/GridLMM_haplotypes/models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',chr,h,pheno,env))
  #heffects=hresults[hresults$X_ID==line$highest_SNP,c(2,6:(h+5))]
  #names(heffects)=c("X_ID",seq(1,h))
  #hmelt=melt(heffects,"X_ID")
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==line$highest_SNP,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$method="Haplotype_probs"
  fmelt$hapgrp=ibd_seg
  names(fmelt)=c('value','se','tvalue','variable','method','hapgrp')

  #names(hmelt)=c('X_ID','hapgrp','value','variable')
  #hmelt$hapgrp=factor(hmelt$hapgrp,order(seq(1,h)))
  #hmelt$method="Haplotype_probs"
  #hmelt=hmelt[,c('X_ID','variable','value','method','hapgrp')]
  fmelt$hapgrp=factor(fmelt$hapgrp,levels=seq(1,h))
  #names(hmelt)=c('value','se','tvalue','variable','method','hapgrp')

  fmelt$variable_f=factor(fmelt$variable,levels=fmelt[order(fmelt$hapgrp),]$variable)

  count=count+1

  h_notf_plots[[count]]=ggplot(fmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge()) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Founder effect Size of %s (QTL in Haplotype notFounder)',q)) + theme(axis.text.x=element_text(size=5),title=element_text(size=10))
  #f_nots_plots[[count]]=ggplot(allmelt,aes(x=variable,y=value)) + geom_bar(aes(fill=allele),stat="identity")+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Highest SNP effect Size of %s (QTL in Founder not SNP)',q))

}

pdf('GridLMM/effect_sizes/Haplotype_notFounder_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(h_notf_plots)){
  print(h_notf_plots[[i]])
}
dev.off()
