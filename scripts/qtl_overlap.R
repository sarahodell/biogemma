#!/usr/bin/env Rscript

library('data.table')
library('tidyr')
library('dplyr')
library('reshape2')


qtl=fread('GridLMM/Biogemma_QTL_all_thresholds.csv',data.table=F)
qtl=qtl[qtl$`10P_Sig`!='T',]
qtl=qtl[qtl$Environment=="ALL",]
#600K SNP
#sqtl=qtl[qtl$Method=="600K_SNP",]
ses=readRDS('GridLMM/effect_sizes/All_QTL_ES.rds')
#n=which(unlist(lapply(s_f_h_data,function(x) x$qtl_id == name & x$focal=="Founder_probs")))

chroms=unique(qtl$Chromosome)
for(c in chroms){
  sub=qtl[qtl$Chromosome==c,]
  qtl_ids=unique(sub$ID)
  ql=length(qtl_ids)
  if(ql>1){
    for(i in 1:(ql-1)){
      for(j in (i+1):ql){
        q1=qtl_ids[i]
        q2=qtl_ids[j]
        env1=as.data.table(qtl[qtl$ID==q1,])
        env2=as.data.table(qtl[qtl$ID==q2,])
        setkey(env2,left_bound_bp,right_bound_bp)
        comparison=foverlaps(env1,env2,by.x=c('left_bound_bp','right_bound_bp'),by.y=c('left_bound_bp','right_bound_bp'),nomatch=NA)
        comparison=comparison[!is.na(comparison$ID),]
        if(dim(comparison)[1]!=0){
          print(sprintf("Overlap between %s and %s",q1,q2))
        }
      }
    }
  }
}

module1=c('qDTS3_2','qHGM3_2','qDTA3_2')
sub3=qtl[qtl$ID %in% module1,]

sub3s=sub3[sub3$Method=="600K_SNP",]

ids=paste0(unique(sub3s$Phenotype),'_ALL_',module1)
dts3_2=ses[[which(unlist(lapply(ses,function(x) x$qtl_id==ids[1] & x$focal=="600K_SNP"& x$method=="600K_SNP")))]]$SE
rownames(dts3_2)=c(1,2)
hgm3_2=ses[[which(unlist(lapply(ses,function(x) x$qtl_id==ids[2] & x$focal=="600K_SNP"& x$method=="600K_SNP")))]]$SE
rownames(hgm3_2)=c(1,2)

dta3_2=ses[[which(unlist(lapply(ses,function(x) x$qtl_id==ids[3] & x$focal=="600K_SNP" & x$method=="600K_SNP")))]]$SE
rownames(dta3_2)=c(1,2)

geno=fread('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr3_binary.csv',data.table=F)
dts_snp=sub3s[sub3s$ID=="qDTS3_2",]$highest_SNP
dta_snp=sub3s[sub3s$ID=="qDTA3_2",]$highest_SNP
hgm_snp=sub3s[sub3s$ID=="qHGM3_2",]$highest_SNP

geno=geno[,c(dts_snp,dta_snp,hgm_snp)]

dts_pop=dts3_2$value[geno[,dts_snp]+1]
dta_pop=dta3_2$value[geno[,dta_snp]+1]
hgm_pop=hgm3_2$value[geno[,hgm_snp]+1]

cor(dts_pop,dta_pop)
#[1] 0.4922273
cor(dts_pop,hgm_pop)
#[1] 0.4029918
cor(dta_pop,hgm_pop)
#[1] 0.5603087

ses=readRDS('GridLMM/effect_sizes/All_effect_sizes.rds')

sub3f=sub3[sub3$Method=="Founder_probs",]
ids=paste0(unique(sub3s$Phenotype),'_ALL_',module1)
dts3_2=ses[[which(unlist(lapply(ses,function(x) x$id==ids[1])))]]$values
hgm3_2=ses[[which(unlist(lapply(ses,function(x) x$id==ids[2])))]]$values
dta3_2=ses[[which(unlist(lapply(ses,function(x) x$id==ids[3])))]]$values
geno=readRDS('genotypes/probabilities/geno_probs/bg3_filtered_genotype_probs.rds')
dts_snp=sub3f[sub3f$ID=="qDTS3_2",]$highest_SNP
dta_snp=sub3f[sub3f$ID=="qDTA3_2",]$highest_SNP
hgm_snp=sub3f[sub3f$ID=="qHGM3_2",]$highest_SNP

geno=lapply(geno,function(x) x[,c(dts_snp,dta_snp,hgm_snp)])
X = do.call(cbind,lapply(geno,function(x) x[,dts_snp]))
dts_pop=X %*% dts3_2$f_value
X = do.call(cbind,lapply(geno,function(x) x[,dta_snp]))
dta_pop=X %*% dta3_2$f_value
X = do.call(cbind,lapply(geno,function(x) x[,hgm_snp]))
hgm_pop=X %*% hgm3_2$f_value


cor(dts_pop,dta_pop)
#[1] 0.9744875
cor(dts_pop,hgm_pop)
#[1] 0.4429077
cor(dta_pop,hgm_pop)
#[1] 0.5603087
