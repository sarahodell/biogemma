library('data.table')
library('dplyr')

X <- fread("out.012", header = FALSE,data.table=F)
X$V1 <- NULL
X <- data.matrix(X)
map <- fread("out.012.pos", header = FALSE,data.table=F)
names(map)=c("chr","pos")
map$snp <- paste0('S',map$chr,'_',map$pos)

Geno <- fread("out.012.indv",header = FALSE,data.table=F)
colnames(X) <- map$snp
rownames(X) <- Geno$V1

phenotypes=fread('phenotypes.csv',data.table=F)
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes <- phenotypes[phenotypes$Genotype_code %in% rownames(X),]

X = X[as.character(phenotypes$Genotype_code),]
X <- unique(X)

# Leave on chromosome out
for(i in 1:10){
  chr=paste0('S',i,'_')
  sub_X = X[,!(grepl(chr,colnames(X)))]
  X_centered<-sweep(sub_X,2,colMeans(sub_X),'-') # center marker genotypes
  K = tcrossprod(X_centered)/ncol(X_centered)
  fwrite(as.data.frame(K),sprintf('K_matrix_chr%.0f.txt',i),row.names=T,quote=F,sep='\t')

}