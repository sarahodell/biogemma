#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')

# Read in Kinship Matrix
K=fread(sprintf('../K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('../phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

# Read in the haplotype group probabilities
# Filter genotypes that are not in the K matrix
X_list=readRDS(sprintf('../../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

mite_prob=fread('../mite_probabilities.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#new_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#Make B73 the first in the list so that it is the one that is dropped
#names(X_list)=founders
#X_list=X_list[new_founders]
# Grab only lines that have the MITE
# Run GridLMM
rownames(mite_prob)=mite_prob$ID
mite_prob=mite_prob[data_blup$ID,]
#rownames(mite_prob)=seq(1,dim(mite_prob)[1])



X_list_order_1=lapply(X_list,function(x) x[data_blup$ID,])
has_mite=mite_prob[mite_prob$final>=0.9,]$ID
X_list_ordered=lapply(X_list_order_1,function(x) x[has_mite,])
data_blup=data_blup[data_blup$ID %in% has_mite,]
K=K[rownames(K) %in% has_mite,colnames(K) %in% has_mite]

Y=as.matrix(data_blup$y)
null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F,verbose=F)

h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup
X_cov=null_model$lmod$X
X_list_null=NULL
gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

saveRDS(gwas,sprintf('models/Biogemma_chr%s_%s_x_ALL_founderprobs_MITE_only.rds',chr,pheno))

# Convert all very high and very low probabilities to 1 and 0, respectively
#X_list_full = lapply(X_list_ordered,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.95,1,ifelse(x[,i]<=0.05,0,x[,i]))))
#for(i in 1:16){dimnames(X_list_full[[i]])[[2]]=dimnames(X_list_ordered[[i]])[[2]]}

#gwas_adjusted=gwas
#sums=lapply(X_list_full,function(x) colSums(x))
#for(i in 1:16){
#    s=sums[[i]]
#    t=dim(X_list_full[[i]])[1]-2
#    l=2
#    grab=which(s>t,s)
#    grab=c(grab,which(s<l,s))
#    grab=sort(grab)
#    beta=sprintf('beta.%.0f',seq(1,16))
#    gwas_adjusted[grab,beta]=0
#    gwas_adjusted[grab,'p_value_ML']=0.99
#    print(grab)
#}

#saveRDS(gwas_adjusted,sprintf('models/Biogemma_chr%s_%s_x_ALL_founderprobs_MITE_only_adjusted.rds',chr,pheno))
