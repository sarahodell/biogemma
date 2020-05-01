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

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite=c(T,F,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

bg=readRDS('genotypes/probabilities/geno_probs/raw/bg8_genoprobs_010319.rds')
pmap=fread('genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)



mite_start=135.947816*1e6
mite_end=135.946644*1e6
region=pmap[pmap$pos>135500000 & pmap$pos<136500000,]$marker
markers=dimnames(bg[[1]])[[3]]
markers=markers[markers %in% region]
sub8=lapply(bg,function(x) x[,,markers])

#Probability that each individual has the MITE
mite_prob=sapply(seq(1,length(region)), function(x) sapply(seq(1,344), function(i) sum(sub8[[1]][i,has_mite,x])))
rownames(mite_prob)=dimnames(bg[[1]])[[1]]
colnames(mite_prob)=region

marker='AX-91102970'

at_mite=as.data.frame(mite_prob[,marker])
at_mite$ID=rownames(at_mite)
names(at_mite)=c(marker,'ID')
at_mite=at_mite[,c('ID',marker)]

fwrite(at_mite,'GridLMM/mite_probabilities.txt',quote=F,row.names=F,sep='\t')
