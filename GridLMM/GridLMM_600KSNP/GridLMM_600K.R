#!/usr/bin/env Rscript

## Run GridLMM using the 600K SNP Array data

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])

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
phenotypes=fread('../phenotypes.csv',data.table=F)
print(dim(phenotypes))
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
print(dim(phenotypes))
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)])
data=data[data$Loc.Year.Treat==env,]

map=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

geno=fread(sprintf('DH_geno_chr%s_121718.csv',chr),data.table=F)
X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
X=t(X)
print(dim(X))
X=as.data.frame(X)
rownames(X)=geno$ind
colnames(X)=colnames(geno)[2:dim(geno)[2]]
X=as.matrix(X)

data=data[data$ID %in% rownames(X),]
#X_list_full=X[data$ID,]



gwas = GridLMM_GWAS(
                        formula = y~1 + (1|ID),
                        test_formula = ~1,
                        reduced_formula = ~0,
                        data = data, 
                        weights = NULL, 
                        X = X,
                        X_ID = 'ID',
                        h2_start = NULL,
                        h2_step = 0.01,
                        max_steps = 10,
                        X_map = map,
                        relmat = list(ID=K),
                        centerX = TRUE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'REML',
                        mc.cores = 1,
                        verbose = FALSE
)

saveRDS(gwas,sprintf('chr%s_%s_x_%s_600KSNP.rds',chr,pheno,env))