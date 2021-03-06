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

print("Reading in phenotypes")
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('../phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]


data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

print("Calculating BLUPs")
m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

#map=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

print("Reading in genotype data")
X=fread(sprintf('../../genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F,stringsAsFactors=F)
#X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
#X=t(X)
#print(dim(X))
#X=as.data.frame(X)
#rownames(X)=geno$ind
#colnames(X)=colnames(geno)[2:dim(geno)[2]]
rownames(X)=X$ind

X=X[,2:dim(X)[2]]
i=intersect(data_blup$ID,rownames(X))
X=X[i,]
data_blup=data_blup[i,]
K=K[i,i]
X=as.matrix(X)
#dimr=dim(X)[1]
#dimc=dim(X)[2]
#mono=c()
#m=apply(X,MARGIN=2,FUN=function(n) length(unique(n)))
#if(length(m[m==dimr]!=0)){
#  mono=c(mono,which(m==dimr))
#}
#if(length(mono)>=1){
#  X_filtered=X[,-mono]
#}else{
#  X_filtered=X
#}

#X=as.matrix(X_filtered)
print("Running GWAS")

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
                        centerX = FALSE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'REML',
                        mc.cores = cores,
                        verbose = FALSE
)


saveRDS(gwas,sprintf('models/chr%s_%s_x_ALL_600KSNP_ML.rds',chr,pheno))

h2=gwas$setup$h2_start
hinfo=data.frame(method="600K_SNP",phenotype=pheno,environment="ALL",chr=chr,h2=h2,hap=NA,stringsAsFactors=F)
fwrite(hinfo,'../heritabilities.txt',quote=F,sep='\t',row.names=F,append=T)
