#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
grp=as.character(args[[4]])
cores=as.numeric(args[[5]])
reps=as.numeric(args[[6]])

library('GridLMM')
library('data.table')
library('dplyr')
library('parallel')
library('MASS')

# Read in Kinship Matrix
K=fread(sprintf('../../K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('../../phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)])
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y = data$y-mean(data$y)

X=fread(sprintf('../../../genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
#X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
#X=t(X)
#print(dim(X))
#X=as.data.frame(X)
rownames(X)=X$ind
#colnames(X)=colnames(geno)[2:dim(geno)[2]]
X=X[,2:dim(X)[2]]
# Drop non-matching IDs
X=X[rownames(X) %in% data$ID,]
data=data[data$ID %in% rownames(X),]

all_gwas=c()

n_reps=seq(1,reps)

randomized_gwas<-function(rep){
   len=dim(X)[1]
   draw=sample(len,len,replace=F)

   X_list_reordered=X[draw,]
   rownames(X_list_reordered)=rownames(X)
   gwas = GridLMM_GWAS(
                        formula = y~1 + (1|ID),
                        test_formula = ~1,
                        reduced_formula = ~0,
                        data = data,
                        weights = NULL,
                        X = X_list_reordered,
                        X_ID = 'ID',
                        h2_start = NULL,
                        h2_step = 0.01,
                        max_steps = 100,
                        relmat = list(ID=K),
                        centerX = TRUE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'ML',
                        mc.cores = 1,
                        verbose = FALSE
			)
   gwas=gwas$results
   gwas=gwas[!is.na(gwas$p_value_ML),]
   tmp=data.frame(chr=chr,rep=rep,pvalue=min(gwas$p_value_ML))
}
print(system.time({
   results=lapply(n_reps,randomized_gwas)
}))

saveRDS(results,sprintf('test_models/chr%s_%s_x_%s_600KSNP_%srep.rds',chr,pheno,env,grp))
