#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')

geno=fread(sprintf('DH_geno_chr%s_121718.csv',chr),data.table=F)
X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
X=t(X)
print(dim(X))
X=as.data.frame(X)
rownames(X)=geno$ind
colnames(X)=colnames(geno)[2:dim(geno)[2]]
X=as.matrix(X)

m_names=colnames(X)
cutoff=0.95
dropped=list()
count=1
dropped[[count]] = list(marker=m_names[1],linked=c())

keep=c(1)
start=unlist(unname(X[,1]))
size=dim(X)[2]
for(i in 2:size){
  ind_max=unlist(unname(X[,i]))
  if(cor(start,ind_max)<cutoff){
    keep=c(keep,i)
    start=ind_max
    start_ind=i
    count=count+1
    dropped[[count]]=list(marker=m_names[start_ind],linked=c())
  }
  else{
    dropped[[count]]$linked=c(dropped[[count]]$linked,m_names[i])
  }
}

X_filtered=X[,keep]

saveRDS(dropped,sprintf('bg%s_dropped_600K.rds',chr))
saveRDS(X_filtered,sprintf('bg%s_filtered_600Kgeno.rds',chr))
sprintf("Keeping %0.f markers on chromosome %s with correlation filter of %.2f",length(keep),chr,cutoff)

