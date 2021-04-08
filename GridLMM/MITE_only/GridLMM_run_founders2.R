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
inds=rownames(X_list[[1]])

mite_prob=fread('../mite_probabilities.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

# Grab only lines that have the MITE
# Run GridLMM
rownames(mite_prob)=mite_prob$ID
mite_prob=mite_prob[data_blup$ID,]
#rownames(mite_prob)=seq(1,dim(mite_prob)[1])
has_mite=mite_prob[mite_prob$final>=0.9,]$ID
i=intercept(has_mite,inds)

X_list_ordered=lapply(X_list,function(x) x[i,])
data_blup=data_blup[i,]
K=K[i,i]

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

#thresholdtable=fread('../threshold_0.05_table.txt',data.table=F)
#cutoff=thresholdtable[thresholdtable$phenotype==pheno & thresholdtable$environment=="ALL" & thresholdtable$method=="founder_probs",]$threshold

#covs=c()
#count=1
#while(max(-log10(gwas$p_value_ML))>=cutoff){
#  cov=gwas[which.max(-log10(gwas$p_value_ML)),]$X_ID
#  covs=c(covs,cov)
#  X_list_ordered=lapply(X_list_ordered,function(x) x[,dimnames(x)[[2]] != cov])
#  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
#  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
#  h2_start
#  V_setup=null_model$setup
#  X_cov=null_model$lmod$X
#  X_list_null=NULL
#  gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
#  count=count+1
#}

#sigs=data.frame(ID=covs,pos=pmap[match(covs,pmap$marker),]$pos,count=seq(1,34),stringsAsFactors=F)

#Use vgt2 state as covariate
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
inds=rownames(X_list[[1]])


zcn8_states=fread('../GridLMM_founderprobs/zcn8_founder_states.txt',data.table=F)
data_blup$vgt2=zcn8_states[match(data_blup$ID,zcn8_states$sample),]$state
data_blup=data_blup[complete.cases(data_blup),]
i=intersect(inds,data_blup$ID)
data_blup=data_blup[i,]
K=K[i,i]
X_list_ordered=lapply(X_list,function(x) x[i,])

Y=as.matrix(data_blup$y)
null_model = GridLMM_ML(y~1 + vgt2 + (1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F,verbose=F)

h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup
X_cov=as.matrix(null_model$lmod$fr[,c('ID','vgt2')])
X_list_null=NULL
gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

saveRDS(gwas,sprintf('models/Biogemma_chr%s_%s_x_ALL_founderprobs_MITE_only.rds',chr,pheno))
