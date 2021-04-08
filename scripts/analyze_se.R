#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('tidyverse')
library('reshape2')
library('lme4')
library('emmeans')
library('multcomp')
library('multcompView')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")

s_f_h_data=readRDS('GridLMM/effect_sizes/All_effect_sizes.rds')
K=fread('GridLMM/K_matrices/K_matrix_chr10.txt',data.table=F)
inds=K[,1]

var_exp=fread('GridLMM/Biogemma_QTL_variance_explained.txt',data.table=F)

qtl_info=c()
for(i in seq(1,length(s_f_h_data))){
  data=s_f_h_data[[i]]$values
  grouptype=s_f_h_data[[i]]$label
  pheno_env_id=s_f_h_data[[i]]$id
  ssnp=s_f_h_data[[i]]$ssnp
  fsnp=s_f_h_data[[i]]$fsnp
  hsnp=s_f_h_data[[i]]$hsnp
  chrom=s_f_h_data[[i]]$chrom

  geno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chrom),data.table=F)
  geno=geno[geno$ind %in% inds,ssnp,drop=F]
  fgeno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%s.csv',chrom),data.table=F)
  smaf=sum(geno[,1])/nrow(geno)
  smaf=c(1-smaf,smaf)
  f_smaf=sum(ifelse(fgeno[,ssnp]=="A",0,1))/16
  f_smaf=c(1-f_smaf,f_smaf)
  if(is.factor(data$allele[1])){
    smaf=smaf[as.numeric(data$allele)]
    f_smaf=f_smaf[as.numeric(data$allele)]
  }
  else{
    smaf=smaf[as.numeric(data$allele) + 1]
    f_smaf=f_smaf[(as.numeric(data$allele) + 1)]
  }
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chrom))
  freps=round(colSums(fprobs[[1]][inds,,fsnp]))
  hreps=round(colSums(fprobs[[1]][inds,,hsnp]))
  names(hreps)=founders
  names(freps)=founders
  hreps2=c()
  fsum=sum(freps)
  hsum=sum(hreps)
  for(h in 1:max(data$hapgrp)){hreps2=c(hreps2,sum(hreps[rownames(data[data$hapgrp==h,])]))}
  hreps=hreps2[data$hapgrp]
  fmaf=freps/fsum
  hmaf=hreps/hsum
  #data
  data$fmaf=fmaf
  data$smaf=smaf
  data$hmaf=hmaf
  data$f_smaf=f_smaf
  pves=var_exp[var_exp$pheno_env_id==pheno_env_id,]
  data$s_pve=pves[pves$method=="600K_SNP",]$pve
  data$f_pve=pves[pves$method=="Founder_probs",]$pve
  data$h_pve=pves[pves$method=="Haplotype_probs",]$pve

  data$svalue_scaled=data$s_value/sqrt(pves[pves$method=="600K_SNP",]$tot_pheno_var)
  data$fvalue_scaled=data$f_value/sqrt(pves[pves$method=="Founder_probs",]$tot_pheno_var)
  data$hvalue_scaled=data$h_value/sqrt(pves[pves$method=="Haplotype_probs",]$tot_pheno_var)
  data$s_se_scaled=data$s_se/sqrt(pves[pves$method=="600K_SNP",]$tot_pheno_var)
  data$f_se_scaled=data$f_se/sqrt(pves[pves$method=="Founder_probs",]$tot_pheno_var)
  data$h_se_scaled=data$h_se/sqrt(pves[pves$method=="Haplotype_probs",]$tot_pheno_var)
  data$grouptype=grouptype
  data$chrom=chrom
  data$pheno_env_id=pheno_env_id
  data$ssnp=ssnp
  data$fsnp=fsnp
  data$hsnp=hsnp
  qtl_info=rbind(qtl_info,data)
}
qtl_info=as.data.frame(qtl_info,stringsAsFactors=F)
rownames(qtl_info)=seq(1,nrow(qtl_info))
#fwrite(qtl_info,'GridLMM/effect_sizes/all_qtl_info.txt',quote=F,sep'\t',row.names=F)
qtl_info=fread('GridLMM/effect_sizes/all_qtl_info.txt',data.table=F)

types=unique(qtl_info$grouptype)

#test1=qtl_info %>% group_by(grouptype) %>% summarize()
h2s=fread('GridLMM/heritabilities.txt',data.table=F)
h2s=h2s[complete.cases(h2s),]
h2s$pheno_env=paste0(h2s$phenotype,'_',h2s$environment)
h2_sum=h2s %>% group_by(pheno_env,method) %>% summarize(h2_avg=mean(ID),h2_sd=sd(ID))

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl_info$phenotype=qtl[match(qtl_info$pheno_env_id,qtl$pheno_env_id),]$Phenotype
qtl_info$environment=qtl[match(qtl_info$pheno_env_id,qtl$pheno_env_id),]$Environment

s_h2s=h2_sum[h2_sum$method=="600K_SNP",]


#qtl_info$phenotype=sapply(seq(1,nrow(qtl_info)),function(x) strsplit(qtl_info))
# For S_only QTL, how do the pves and std. errors compare to S_F_H
#test1=qtl_info[qtl_info$grouptype =="S_F_H" | qtl_info$grouptype == "S_only",]
t1=unique(qtl_info$pheno_env_id)
a=c()
for(t in t1){
  data=qtl_info[qtl_info$pheno_env_id==t,]
  type=unique(data$grouptype)
  svalue0=as.numeric(unique(data[data$allele==0,]$svalue_scaled))
  svalue1=as.numeric(unique(data[data$allele==1,]$svalue_scaled))
  sse0=as.numeric(unique(data[data$allele==0,]$s_se_scaled))
  sse1=as.numeric(unique(data[data$allele==1,]$s_se_scaled))
  spve=as.numeric(unique(data[data$allele==1,]$s_pve))
  smaf0=as.numeric(unique(data[data$allele==0,]$smaf))
  if(smaf0>0.5){
    smaf0=1-smaf0
  }
  #print(smaf0)
  f_smaf=as.numeric(unique(data[data$allele==0,]$f_smaf))
  if(f_smaf>0.5){
    f_smaf=1-f_smaf
  }
  a=rbind(a,c(t,type,svalue0,svalue1,sse0,sse1,spve,smaf0,f_smaf))
}
a=as.data.frame(a,,stringsAsFactors=F)
names(a)=c('pheno_env_id','grouptype','svalue0','svalue1','sse0','sse1','spve','smaf0','f_smaf')
a[, c(3:9)] <- sapply(a[, c(3:9)], as.numeric)

a$value_diff=abs(a$svalue0-a$svalue1)
a$avg_se = apply(a[,c('sse1','sse0')],MARGIN=1,mean)
in_s=c('S_only','S_F_H','S_and_F')
a$snp_found = ifelse(a$grouptype %in% in_s,T,F)
#m1=lmer(value_diff ~ (1|grouptype) ,a)
m1=lm(value_diff ~ 0 + snp_found,a)
es_emmeans=emmeans(m1,'snp_found')
comp1=contrast(es_emmeans)
#No

m2 = lm(avg_se ~ 0 + snp_found,a)
se_emmeans=emmeans(m2,'snp_found')
comp2=contrast(se_emmeans)
#Standard Errors are larger in S_only than in those found in S_F_H
#summary(m2)$coef
#Yes!
#Analysis of Variance Table

#Response: avg_se
#          Df   Sum Sq  Mean Sq F value    Pr(>F)
#grouptype  2 0.281166 0.140583  922.59 < 2.2e-16 ***
#Residuals 30 0.004571 0.000152

m3 = lm(smaf0 ~ 0 + snp_found,a)
maf_emmeans=emmeans(m3,'snp_found')
comp3=contrast(maf_emmeans)
# Slightly significant (p = 0.0479)
# Higher maf in thos found by S than not

#summary(m3)$coef

m4 = lm(spve ~ 0 + snp_found,a)
pve_emmeans=emmeans(m4,'snp_found')
comp4=contrast(pve_emmeans)
# not significantly different pve for those found by snps and those not

m5 = lm(f_smaf ~ 0 + snp_found,a)
f_smaf_emmeans=emmeans(m5,'snp_found')
comp5=contrast(f_smaf_emmeans)
# F only vs S_F_H

a=c()
for(t in t1){
  data=qtl_info[qtl_info$pheno_env_id==t,]
  type=unique(data$grouptype)
  fvalues=c()
  for(f in founders){
    fvalue=as.numeric(unique(data[data$founder==f,]$fvalue_scaled))
    fse=as.numeric(unique(data[data$founder==f,]$f_se_scaled))
    fpve=as.numeric(unique(data[data$founder==f,]$f_pve))
    #sse1=as.numeric(unique(data[data$allele=="1",]$s_se_scaled))
    fmaf=as.numeric(unique(data[data$founder==f,]$fmaf))
    a=rbind(a,c(t,f,type,fvalue,fse,fpve,fmaf))
  }
}
a=as.data.frame(a,,stringsAsFactors=F)
names(a)=c('pheno_env_id','founder','grouptype','fvalue','fse','fpve','fmaf')
a[, c(4:7)] <- sapply(a[, c(4:7)], as.numeric)

#a$value_diff=abs(a$svalue0-a$svalue1)
#a$avg_se = apply(a[,c('sse1','sse0')],MARGIN=1,mean)
in_f=c('F_only','S_F_H','F_and_H')
a$f_found = ifelse(a$grouptype %in% in_f,T,F)

#m1=lmer(value_diff ~ (1|grouptype) ,a)
#m1=lm(value_diff ~ 0 + grouptype,a)
#es_emmeans=emmeans(m1,'grouptype')
#comp1=contrast(es_emmeans)
#No
m1 = lm(fmaf ~ 0 + f_found,a)
fmaf_emmeans=emmeans(m1,'f_found')
comp3=contrast(maf_emmeans)


m2 = lm(fse ~ 0 + f_found,a)
se_emmeans=emmeans(m2,'f_found')
comp2=contrast(se_emmeans)

m3 = lm(fmaf ~ 0 + f_found,a)
maf_emmeans=emmeans(m3,'f_found')
comp3=contrast(maf_emmeans)

m4 = lm(fpve ~ 0 + f_found,a)
pve_emmeans=emmeans(m4,'f_found')
comp4=contrast(pve_emmeans)

 #significant differenc in pve***
 # higher pve for those found by F than not
#summary(m3)$coef

# H only

a=c()
for(t in t1){
  data=qtl_info[qtl_info$pheno_env_id==t,]
  type=unique(data$grouptype)
  hvalues=c()
  maxh=max(data$hapgrp)
  for(h in 1:maxh){
    hvalue=as.numeric(unique(data[data$hapgrp==h,]$hvalue_scaled))
    hse=as.numeric(unique(data[data$hapgrp==h,]$h_se_scaled))
    hpve=as.numeric(unique(data[data$hapgrp==h,]$h_pve))
    #sse1=as.numeric(unique(data[data$allele=="1",]$s_se_scaled))
    hmaf=as.numeric(unique(data[data$hapgrp==h,]$hmaf))
    if(hmaf>0.5){
      hmaf=1-hmaf
    }
    a=rbind(a,c(t,h,type,hvalue,hse,hpve,hmaf))
  }
}
a=as.data.frame(a,,stringsAsFactors=F)
names(a)=c('pheno_env_id','hapgrp','grouptype','hvalue','hse','hpve','hmaf')
a[, c(2,4:7)] <- sapply(a[, c(2,4:7)], as.numeric)

#a$value_diff=abs(a$svalue0-a$svalue1)
#a$avg_se = apply(a[,c('sse1','sse0')],MARGIN=1,mean)
in_h=c('H_only','S_F_H','F_and_H')
a$h_found = ifelse(a$grouptype %in% in_h,T,F)

#m1=lmer(value_diff ~ (1|grouptype) ,a)
#m1=lm(value_diff ~ 0 + grouptype,a)
#es_emmeans=emmeans(m1,'grouptype')
#comp1=contrast(es_emmeans)
#No
a2 = a %>% group_by(pheno_env_id,h_found) %>% summarize(hmax=max(hapgrp))
m1=lm(hmax ~ 0 + h_found,a2)
hmax_emmeans=emmeans(m1,'h_found')
comp1=contrast(hmax_emmeans)

m2 = lm(hse ~ 0 + h_found,a)
se_emmeans=emmeans(m2,'h_found')
comp2=contrast(se_emmeans)

m3 = lm(hmaf ~ 0 + h_found,a)
maf_emmeans=emmeans(m3,'h_found')
comp3=contrast(maf_emmeans)

m4 = lm(hpve ~ 0 + h_found,a)
pve_emmeans=emmeans(m4,'h_found')
comp4=contrast(pve_emmeans)
# pve is significant <0.0001
# higher pve explained for the ones that H found

var_exp$contrast=paste0(var_exp$method,var_exp$sig)
m6=lmer(pve ~ 0 + contrast + (1|pheno_env_id), var_exp)
emmeans6=emmeans(m6,'contrast')
comp6=as.data.frame(contrast(emmeans6,'pairwise'),stringsAsFactors=)
emmeans6=as.data.frame(cld(emmeans6,Letters = letters),stringsAsFactors=F)
emmeans6=emmeans6[order(emmeans6$.group),]
emmeans6$found=c(T,T,F,F,F,T)

png('pve_contrasts.png',width=800,height=800)
print(ggplot(emmeans6,aes(x=contrast,y=emmean,color=found)) + geom_point() +
geom_text(aes(label=.group),hjust=2,vjust=0) +
geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) +
xlab("Comparisons") +
ylab("Percent Phenotypic Variance Explained"))
dev.off()
#contrast                                   estimate      SE    df t.ratio p.value
# 600K_SNPFALSE - 600K_SNPTRUE                0.02337 0.00715 107.7   3.267  0.0178
# 600K_SNPFALSE - Founder_probsFALSE         -0.08007 0.00912 107.3  -8.782  <.0001
# 600K_SNPFALSE - Founder_probsTRUE          -0.11577 0.00677 100.7 -17.112  <.0001
# 600K_SNPFALSE - Haplotype_probsFALSE       -0.07920 0.00839 105.5  -9.446  <.0001
# 600K_SNPFALSE - Haplotype_probsTRUE        -0.10852 0.00686 101.3 -15.822  <.0001
# 600K_SNPTRUE - Founder_probsFALSE          -0.10345 0.00750 106.0 -13.792  <.0001
# 600K_SNPTRUE - Founder_probsTRUE           -0.13915 0.00445  92.3 -31.278  <.0001
# 600K_SNPTRUE - Haplotype_probsFALSE        -0.10258 0.00665 104.2 -15.425  <.0001
# 600K_SNPTRUE - Haplotype_probsTRUE         -0.13189 0.00457  93.2 -28.883  <.0001
# Founder_probsFALSE - Founder_probsTRUE     -0.03570 0.00782 113.5  -4.566  0.0002
# Founder_probsFALSE - Haplotype_probsFALSE   0.00087 0.00826  91.8   0.105  1.0000
# Founder_probsFALSE - Haplotype_probsTRUE   -0.02845 0.00783 113.3  -3.632  0.0055
# Founder_probsTRUE - Haplotype_probsFALSE    0.03657 0.00683 110.8   5.356  <.0001
# Founder_probsTRUE - Haplotype_probsTRUE     0.00725 0.00438  87.2   1.658  0.5629
# Haplotype_probsFALSE - Haplotype_probsTRUE -0.02932 0.00707 113.5  -4.145  0.0009

#Coefficient of Variation (SE/ES)
tmp=qtl_info %>% group_by(pheno_env_id) %>% summarize(h_cv=mean(h_se_scaled),f_cv=mean(f_se_scaled),s_cv=mean(s_se_scaled))
tmelt=melt(tmp,'pheno_env_id')
contrast=c()
for(i in 1:nrow(tmelt)){
  line=tmelt[i,]
  name=line$pheno_env_id
  method=line$variable
  if(method=="h_cv"){
    m="Haplotype_probs"
  }
  else if(method=="f_cv"){
    m="Founder_probs"
  }
  else{
    m="600K_SNP"
  }
  sig=var_exp[var_exp$pheno_env_id==name & var_exp$method==m,]$sig
  contrast=c(contrast,paste0(m,sig))
}
tmelt$contrast=contrast
m8=lmer(value ~ 0 + contrast + (1|pheno_env_id), tmelt)
emmeans8=emmeans(m8,'contrast')
comp8=as.data.frame(contrast(emmeans8,'pairwise'),stringsAsFactors=)
emmeans8=as.data.frame(cld(emmeans8,Letters = letters),stringsAsFactors=F)
emmeans8=emmeans8[order(emmeans8$.group),]
emmeans6$found=c(T,T,F,F,F,T)
