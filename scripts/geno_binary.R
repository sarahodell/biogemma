#/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')
library('tibble')

geno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_121718.csv',chr),data.table=F)
X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
X=t(X)
X=as.data.frame(X,stringsAsFactors=F)
rownames(X)=geno$ind
colnames(X)=colnames(geno)[2:dim(geno)[2]]
X<-rownames_to_column(X,"ind")
X=X[,c("ind",colnames(geno)[2:dim(geno)[2]])]
fwrite(X,sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),quote=F,sep=',')
