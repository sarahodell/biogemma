#!/usr/bin/env Rscript

library('data.table')
library('rrBLUP')
library('lme4')
#all_genos=c()
#for(chr in 1:10){
#  geno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%.0f_121718.csv',chr),data.table=F)
#  geno=t(geno)
#  geno=as.data.frame(geno,stringsAsFactors=F)
#  inds=unlist(unname(geno[1,]))
#  markers=unlist(unname(rownames(geno)[-1]))
#  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
#  geno=geno[-1,]
#  geno=apply(geno,MARGIN=2,function(x) ifelse(x=="A",-1,1))
#  geno=as.data.frame(geno,stringsAsFactors=F)
#  names(geno)=inds
#  geno$marker=markers

#  geno$chr=chr
#  geno$pos=pmap[match(geno$marker,pmap$marker),]$pos
#  geno=geno[,c('marker','chr','pos',inds)]
#  geno=geno[order(geno$pos),]
#  all_genos=rbind(all_genos,geno)
#}

#all_genos=as.data.frame(all_genos,stringsAsFactors=F)
#names(all_genos)=c('marker','chr','pos',inds)
#fwrite(all_genos,'genotypes/qtl2/Biogemma_DHgenos/rrBLUP_genotype.txt',row.names=F,quote=F,sep=',')

pheno="female_flowering_d6"
phenotypes=fread('GridLMM/phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

print("Calculating BLUPs")
m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

keptsnps=fread('genotypes/600K/Biogemma_LD_0.2_filtered_SNPs.txt',data.table=F,header=F)
genos=fread('genotypes/qtl2/Biogemma_DHgenos/rrBLUP_genotype.txt',data.table=F)


#genos=genos[genos$marker %in% keptsnps$V1,]
allres=c()
#gv_list=list()
us=c()
for(chr in 1:10){
  K=fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-" ,".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  geno=genos[genos$chr==chr,]

  i=intersect(data_blup$ID,names(geno)[4:347])
  geno=geno[,c('marker','chr','pos',i)]
  data_blup=data_blup[data_blup$ID %in% i,]
  K=K[i,i]
## TO DO: SNP Pruning with snpgdsLDpruning {SNPRelate}
  X=geno[,i]
  colnames(X)=i
  rownames(X)=geno$marker
  X=t(X)
  #index=names(which(duplicated(X, MARGIN = 1)))
  #K=K[rownames(K)!=index,colnames(K)!=index]
  #y=data_blup[data_blup$ID!=index,c('y')]
  #X=[rownames(X)!=index,]
  #Zs <- scale(X, center = TRUE, scale = TRUE)
# dimensions
  #n <- nrow(Zs)
  #m <- ncol(Zs)
  #fit1<- mixed.solve(y=data_blup$y,K=K)
  fit2 <- mixed.solve(y = data_blup$y, Z=X)

  u=as.matrix(fit2$u)

  res=X %*% u
  u=as.data.frame(u,stringsAsFactors=F)
  names(u)=c('effect')
  u$chr=chr
  u$marker=rownames(u)
  us=rbind(us,u)
  allres=cbind(allres,res)
  # marker additive genetic variance
  #fit2$Vu
  # residual variance
  #fit2$Ve
  # intercept
  #fit2$beta
  # marker additive genetic effects
  #head(fit2$u)
  #tail(fit2$u)
  # ratio of variance components
  #fit2$Ve / fit2$Vu
  #gv=kin.blup(data=data_blup,geno="ID",pheno='y',GAUSS=FALSE,K=K,fixed=NULL,covariate=NULL,
  #PEV=FALSE,n.core=1,theta.seq=NULL)
  #gv_list[[chr]]=gv$pred
}
pgs=rowSums(allres)
pgs=as.data.frame(pgs,stringsAsFactors=F)
pgs$ind=rownames(pgs)
fwrite(pgs,'GridLMM/pgs/Biogemma_DTS_FT_PGS.txt',row.names=F,quote=F,sep='\t')

us=as.data.frame(us,stringsAsFactors=F)
#names(us)=c('marker','effect','chr')
fwrite(us,'run_magicsim/pgs/DTS_marker_effect_sizes.txt',row.names=F,quote=F,sep='\t')
