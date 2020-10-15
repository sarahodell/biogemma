#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')

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

data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y = data$y - mean(data$y)

# Read in the haplotype group probabilities
# Filter genotypes that are not in the K matrix
X_list=readRDS(sprintf('../../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

new_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite_prob=fread('../mite_probabilities.txt',data.table=F)
rownames(mite_prob)=mite_prob$ID
mite_prob=mite_prob[data$ID,]
rownames(mite_prob)=seq(1,dim(mite_prob)[1])
#Make B73 the first in the list so that it is the one that is dropped
names(X_list)=founders
X_list=X_list[new_founders]

X_list_order_1=lapply(X_list,function(x) x[data$ID,])
has_mite=mite_prob[mite_prob$`AX-91102970`>=0.9,]$ID
X_list_ordered=lapply(X_list_order_1,function(x) x[has_mite,])


size=dim(X_list_ordered[[1]])[2]
#Remove sites with low founder representation
# drop sites with summed founder prob of less than 1
f_sums=lapply(X_list_ordered,function(x) colSums(x))
low_rep=lapply(f_sums,function(x) which(x<1,arr.ind=T))
#low_rep=as.data.frame(low_rep,stringsAsFactors=F)
#which_f=unique(low_rep$row)

dropped=vector("list",length=16)
for(f in 1:16){
  cols=names(low_rep[[f]])
  #geno[[1]][,f,cols]=NA
  dropped[[f]]=cols
}

f_sums2=lapply(X_list_ordered,function(x) colSums(x>=0.8))
low_rep2=lapply(f_sums2,function(x) which(x<5,arr.ind=T))
#ow_rep2=as.data.frame(low_rep2,stringsAsFactors=F)
#which_f2=unique(low_rep2$row)
for(f in 1:16){
  cols=names(low_rep2[[f]])
  #geno[[1]][,f,cols]=NA
  dropped[[f]]=unique(c(dropped[[f]],cols))
}

saveRDS(dropped,sprintf('../../genotypes/probabilities/geno_probs/dropped/bg%s_low_rep_markers_MITE_only.rds',chr))


data=data[data$ID %in% has_mite,]
Y=as.matrix(data$y)

X_list_ordered=lapply(X_list,function(x) x[data$ID,])
# Run GridLMM
null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F)

h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup
X_cov=null_model$lmod$X
X_list_null=NULL

gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

saveRDS(gwas,sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs_MITE_only.rds',chr,pheno,env))
