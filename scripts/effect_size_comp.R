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
qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)

f_ses=readRDS('GridLMM/effect_sizes/Founder_prob_QTL_SEs.rds')
h_ses=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
# Founder vs. Haplotype for shared QTL
q1=qtl_overlap[!is.na(qtl_overlap$S_F_H) | !is.na(qtl_overlap$F_and_H),]
rownames(q1)=seq(1,nrow(q1))

m1=c('Founder_probs','Haplotype_probs')
f_h = qtl[qtl$pheno_env %in% q1$pxe & qtl$Method %in% m1,]
f_h$pheno_env_id = paste0(f_h$pheno_env, '_',f_h$ID)
counts=f_h %>% group_by(pheno_env_id) %>% count
keep=counts[counts$n>=2,]$pheno_env_id
f_h = f_h[f_h$pheno_env_id %in% keep,]
rownames(f_h)=seq(1,nrow(f_h))

f_h_ids=unique(f_h$pheno_env_id)

f_h_plots=list()
count=0

for(q in f_h_ids){
  sub=f_h[f_h$pheno_env_id == q,]
  rownames(sub)=seq(1,nrow(sub))
  for(i in 1:nrow(sub)){
    line=sub[i,]
    pheno=line$Phenotype
    env=line$Environment
    chr=line$Chromosome
    if(line$Method==m1[2]){
      hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
      h=hap_table[hap_table$SNP == line$highest_SNP, ]$HAPGRP
      #hresults=readRDS(sprintf('GridLMM/GridLMM_haplotypes/models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',chr,h,pheno,env))
      #heffects=hresults[hresults$X_ID==line$highest_SNP,c(2,6:(h+5))]
      #names(heffects)=c("X_ID",seq(1,h))
      hmelt=as.data.frame(h_ses[[which(sapply(h_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)

      #hmelt=melt(heffects,"X_ID")
      ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
      pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
      pos=pmap[pmap$marker==line$highest_SNP,]$pos
      ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
      hmelt=hmelt[ibd_seg,]
      rownames(hmelt)=seq(1,nrow(hmelt))
      hmelt$variable=founders
      hmelt$method="Haplotype_probs"
      #names(hmelt)=c('X_ID','hapgrp','value','variable')
      hmelt$hapgrp=factor(ibd_seg,levels=seq(1,16))

      names(hmelt)=c('value','se','tvalue','variable','method','hapgrp')

      hmelt$variable_f=factor(hmelt$variable,levels=hmelt[order(hmelt$hapgrp),]$variable)

    }
    else{
      fmelt=as.data.frame(f_ses[[which(sapply(f_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
      #fresults=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
      #feffects=fresults[fresults$X_ID==line$highest_SNP,c(2,6:21)]
      #names(feffects)=c('X_ID',founders)
      #fmelt=melt(feffects,'X_ID')
      fmelt$variable=founders
      rownames(fmelt)=founders
      fmelt$method="Founder_probs"
      #fmelt$hapgrp=factor(seq(1,16),order(seq(1,16)))
      #fmelt$method="Founder_probs"
      fmelt$hapgrp=factor(seq(1,16),levels=seq(1,16))
      #fmelt$hapgrp=factor(hmelt$hapgrp,levels=seq(1,h))
      names(fmelt)=c('value','se','tvalue','variable','method','hapgrp')
      fmelt$variable_f=factor(fmelt$variable,levels=fmelt[order(fmelt$hapgrp),]$variable)
    }
  }
  allmelt=rbind(hmelt,fmelt)
  count=count+1
  f_h_plots[[count]]=ggplot(allmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))
}

pdf('GridLMM/effect_sizes/Founder_Haplotype_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(f_h_plots)){
  print(f_h_plots[[i]])
}
dev.off()

# Founder vs. SNP overlap for only shared SNPs

q2=qtl_overlap[!is.na(qtl_overlap$S_F_H) | !is.na(qtl_overlap$S_and_F),]
rownames(q2)=seq(1,nrow(q2))

m2=c('Founder_probs','600K_SNP')
f_s = qtl[qtl$pheno_env %in% q2$pxe & qtl$Method %in% m2,]
f_s$pheno_env_id = paste0(f_s$pheno_env, '_',f_s$ID)
counts=f_s %>% group_by(pheno_env_id) %>% count
keep=counts[counts$n>=2,]$pheno_env_id
f_s = f_s[f_s$pheno_env_id %in% keep,]
rownames(f_s)=seq(1,nrow(f_s))

f_s_ids=unique(f_s$pheno_env_id)

f_s_plots=list()
count=0

for(q in f_s_ids){
  sub=f_s[f_s$pheno_env_id == q,]
  rownames(sub)=seq(1,nrow(sub))
  for(i in 1:nrow(sub)){
    line=sub[i,]
    pheno=line$Phenotype
    env=line$Environment
    chr=line$Chromosome
    if(line$Method==m2[2]){
      sresults=readRDS(sprintf('GridLMM/GridLMM_600KSNP/models/chr%.0f_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
      seffects=sresults$results[sresults$results$X_ID==line$highest_SNP,c(2,6:7)]
      names(seffects)=c('X_ID',seq(0,1))
      smelt=melt(seffects,'X_ID')
      smelt$method="600K_SNP"
      ssnp=line$highest_SNP
      geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f_121718.csv',chr),data.table=F)
      rownames(geno)=geno$ind
      alleles=ifelse(geno[founders,ssnp]=="A",1,2)
      smelt=smelt[alleles,]
      rownames(smelt)=seq(1,nrow(smelt))
      smelt$founders=founders
      smelt$variable=as.character(smelt$variable)
      names(smelt)=c('X_ID','allele','value','method','variable')
      smelt=smelt[,c('X_ID','variable','value','method','allele')]
    }
    else{
      fresults=readRDS(sprintf('GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
      feffects=fresults[fresults$X_ID==line$highest_SNP,c(2,6:21)]
      names(feffects)=c('X_ID',founders)
      fmelt=melt(feffects,'X_ID')
      fmelt$method="Founder_probs"
      fmelt$allele=as.character(seq(3:18))
    }
  }
  allmelt=rbind(smelt,fmelt)
  count=count+1
  f_s_plots[[count]]=ggplot(allmelt,aes(x=variable,y=value)) + geom_bar(aes(fill=allele),stat="identity")+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Allele") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr))
}

pdf('GridLMM/effect_sizes/Founder_SNPGWAS_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(f_s_plots)){
  print(f_s_plots[[i]])
}
dev.off()

#Founder vs. SNP overlap for non-shared SNPs
