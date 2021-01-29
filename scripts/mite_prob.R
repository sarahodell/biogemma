#!/usr/bin/env Rscript

pheno="male_flowering_d6"
chr="8"

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('ggplot2')
library('reshape2')
library('tibble')
library('dplyr')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

bg=readRDS('genotypes/probabilities/geno_probs/raw/bg8_genoprobs_010319.rds')
pmap=fread('genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)


founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
names(founder_probs)=founders
mite_start=135.947816*1e6
mite_end=135.946644*1e6
region=pmap[pmap$pos>135500000 & pmap$pos<136500000,]$marker
markers=dimnames(founder_probs[[1]])[[2]]
markers=markers[markers %in% region]
snp="AX-91102858"
sub8=lapply(founder_probs,function(x) x[,markers])
l=dim(sub8[[1]])[1]
mitefounders=founder_probs[c(has_mite)]
mite_prob = rowSums(do.call(cbind,lapply(mitefounders,function(x) x[,markers[1]])))
#mite_prob=ifelse(mite_prob>=0.85,1,0)
mite_prob2 = rowSums(do.call(cbind,lapply(mitefounders,function(x) x[,markers[2]])))

t=data.frame(ID=names(mite_prob),m1=mite_prob,m2=mite_prob2,stringsAsFactors=F)
t$final=rowMeans(t[,c('m1','m2')])

t=t[,c('ID','final')]

#mite_prob2=ifelse(mite_prob2>=0.85,1,0)
#mite_prob=rowSums(lapply(mitefounders,function(x) x))
#Probability that each individual has the MITE
#mite_prob=sapply(seq(1,length(region)), function(x) sapply(seq(1,l), function(i) sum(sub8[i,has_mite,x])))
#rownames(mite_prob)=dimnames(founder_probs[[1]])[[1]]
#colnames(mite_prob)=region

marker='AX-91102970'

at_mite=as.data.frame(mite_prob[,marker])
at_mite$ID=rownames(at_mite)
names(at_mite)=c(marker,'ID')
at_mite=at_mite[,c('ID',marker)]

fwrite(at_mite,'GridLMM/mite_probabilities.txt',quote=F,row.names=F,sep='\t')
