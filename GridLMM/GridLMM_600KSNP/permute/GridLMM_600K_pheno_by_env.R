#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])
reps=as.numeric(args[[5]])

library('GridLMM')
library('data.table')
library('dplyr')


# Read in Kinship Matrix
K=fread(sprintf('../../K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('../../phenotypes.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)])
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y = data$y-mean(data$y)

geno=fread(sprintf('../DH_geno_chr%s_121718.csv',chr),data.table=F)
X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
X=t(X)
print(dim(X))
X=as.data.frame(X)
rownames(X)=geno$ind
colnames(X)=colnames(geno)[2:dim(geno)[2]]

# Drop non-matching IDs
X=X[rownames(X) %in% data$ID,]
data=data[data$ID %in% rownames(X),]
X=as.matrix(X)

all_gwas=c()

for(i in 1:reps){
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
                        method = 'REML',
                        mc.cores = 1,
                        verbose = FALSE
			)
   gwas=gwas$results
   gwas=gwas[!is.na(gwas$p_value_REML),]
   all_gwas=rbind(all_gwas,data.frame(chr=chr,rep=i,pvalue=min(gwas$p_value_REML))) 
}
fwrite(all_gwas,sprintf('test_models/chr%s_%s_x_%s_600KSNP_%.0frep.txt',chr,pheno,env,reps),quote=F,row.names=F)
