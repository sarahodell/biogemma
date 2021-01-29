#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('tidyverse')
library('lme4qtl')

chr="10"


X = fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
rownames(X)=X$ind
x="AX-91157777"
founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

X=X[rownames(X) %in% rownames(founder_probs[[1]]),]
X = X[,x,drop=F]
founder_probs=lapply(founder_probs,function(i) i[,x,drop=F])
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#ref=c("OH43_inra","A632_usa","FV252_inra","F492","B96","CO255_inra","B73_inra","D105_inra")

fgeno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%s_121718.csv',chr),data.table=F)
fgeno=fgeno[,c('ind',x)]
fgeno[,x] = ifelse(fgeno[,x]=="A",0,1)
ref=fgeno[fgeno[,x]==0,]$ind
alt=fgeno[fgeno[,x]==1,]$ind
effects=c()
for(f in founders){effects=c(effects,ifelse(f %in% ref,0,5))}
l=length(founder_probs[[1]])
founder_phenos=sapply(seq(1,l),function(j) sum(unlist(lapply(founder_probs,function(i) i[j,]))*effects))
phenotype=data.frame(ID=rownames(founder_probs[[1]]),y=founder_phenos)
phenotype$y=phenotype$y + rnorm(l)


K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

phenotype = subset(phenotype,ID %in% rownames(K))
fX = do.call(cbind,lapply(founder_probs,function(x) x[,1]))
colnames(fX) = founders
rownames(fX) = dimnames(founder_probs[[1]])[[1]]

#SNPs
X = X[rownames(X) %in% phenotype$ID,]
phenotype$X=X[phenotype$ID]
m0 = relmatLmer(y ~ 0 + X + (1|ID),data=phenotype,relmat = list(ID=K))
s_se=as.data.frame(summary(m0)$coef)

#Founders
fX = fX[rownames(fX) %in% phenotype$ID,]
m1 = relmatLmer(y ~ 0 + fX + (1|ID),data=phenotype,relmat = list(ID=K))
f_se=as.data.frame(summary(m1)$coef)


#Now do with GridLMM...
