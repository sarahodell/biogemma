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

#map=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

X=fread(sprintf('../../../Biogemma_DHgenos/bg%s_filtered_600K.csv',chr),data.table=F,stringsAsFactors=F)
#X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
#X=t(X)
#print(dim(X))
#X=as.data.frame(X)
#rownames(X)=geno$ind
#colnames(X)=colnames(geno)[2:dim(geno)[2]]
rownames(X)=X$ind
X=X[,2:dim(X)[2]]

X=X[rownames(X) %in% data_blup$ID,]
X=as.matrix(X)
data_blup=data_blup[data_blup$ID %in% rownames(X),]

dimr=dim(X)[1]
#dimc=dim(X)[2]
mono=c()
m=apply(X,MARGIN=2,FUN=function(n) length(unique(n)))
if(length(m[m==dimr]!=0)){
  mono=c(mono,which(m==dimr))
}
if(length(mono)>=1){
  X_filtered=X[,-mono]
}else{
  X_filtered=X
}

X=as.matrix(X_filtered)


gwas = GridLMM_GWAS(
                        formula = y~1 + (1|ID),
                        test_formula = ~1,
                        reduced_formula = ~0,
                        data = data_blup,
                        weights = NULL,
                        X = X,
                        X_ID = 'ID',
                        h2_start = NULL,
                        h2_step = 0.01,
                        max_steps = 100,
                        relmat = list(ID=K),
                        centerX = TRUE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'ML',
                        mc.cores = cores,
                        verbose = FALSE
)

saveRDS(gwas,sprintf('models/chr%s_%s_x_ALL_600KSNP_ML.rds',chr,pheno))
