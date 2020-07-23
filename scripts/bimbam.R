#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')

df=fread(cmd=sprintf("gzip -dc bg%s_wgs_alleleprobs.txt.gz",chr),data.table=F)
#df=fread('test.txt',data.table=F)
K=fread(sprintf('../../../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

#Drop samples without phenotype data
df=df[,c('marker','alt1','ref',rownames(K))]

freq=rowSums(df[4:dim(df)[2]])/(dim(df)[2]-3)
min_allele=freq<=0.5

is_major=which(min_allele==F)

bimbam=function(i,df){
  line=df[i,]
  tmp=df[i,'alt1']
  df[i,'alt1']=df[i,'ref']
  df[i,'ref']=tmp
  df[i,4:dim(df)[2]]=1-df[i,4:dim(df)[2]]
  return(df)
}


df_edit <- bimbam(is_major,df)

#for(i in is_major){
#  line=df[i,]
#  tmp=df[i,'alt1']
#  df[i,'alt1']=df[i,'ref']
#  df[i,'ref']=tmp
#  df[i,4:dim(df)[2]]=1-df[i,4:dim(df)[2]]
#}
#t=sapply(is_major,function(x) return(bimbam(x,df)))
df_edit[,4:dim(df_edit)[2]] = df_edit[,4:dim(df_edit)[2]]*2

#fwrite(df,'test_bimbam.txt',row.names=F,quote=F,sep=',',col.names=F)
fwrite(df_edit,sprintf('bg%s_wgs_alleleprobs_bimbam.txt',chr),row.names=F,quote=F,col.names=F,sep=',')
#system(sprintf("gzip bg%s_wgs_alleleprobs_bimbam.txt",chr))
