#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
rep=as.numeric(args[[1]])

#library('SNPRelate')
library('data.table')
library('rrBLUP')
library('lme4')

sample.id=paste0('Sim',seq(1,344))
pop_code=c(rep("DH",344))


#vcf.fn<-sprintf('merged_vcfs/MAGIC_DHSimAll_rep%.0f.vcf.gz',rep)
#snpgdsVCF2GDS(vcf.fn,sprintf("pgs/Rep%.0f.gds",rep),method="biallelic.only")
#genofile<-snpgdsOpen(sprintf("pgs/Rep%.0f.gds",rep))
#snpgdsSummary(sprintf("pgs/Rep%.0f.gds",rep))

#snpset1<-snpgdsLDpruning(genofile,ld.threshold = 0.2)
#With LD 0.2 used 18,423 markers
#snpset1.id<-unlist(snpset1)

#snprslist=read.gdsn(index.gdsn(genofile,"snp.rs.id"))[snpset1.id]


all_genos=c()
for(chr in 1:10){
  geno=fread(sprintf('qtl2_files/MAGIC_DHSimAll_rep%.0f_chr%.0f.csv',rep,chr),data.table=F)
  markers=names(geno)[-1]
#  rownames(geno)=sample.id
#  geno=geno[,snprslist]
  geno=t(geno)
  geno=as.data.frame(geno,stringsAsFactors=F)
  inds=sample.id
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  geno=geno[-1,]
  geno=apply(geno,MARGIN=2,function(x) ifelse(x=="A",-1,1))
  geno=as.data.frame(geno,stringsAsFactors=F)
  names(geno)=inds
  geno$marker=markers

  geno$chr=chr
  geno$pos=pmap[match(geno$marker,pmap$marker),]$pos
  geno=geno[,c('marker','chr','pos',inds)]
  geno=geno[order(geno$pos),]
  all_genos=rbind(all_genos,geno)
}

all_genos=as.data.frame(all_genos,stringsAsFactors=F)
names(all_genos)=c('marker','chr','pos',inds)
#fwrite(all_genos,'genotypes/qtl2/Biogemma_DHgenos/rrBLUP_genotype.txt',row.names=F,quote=F,sep=',')

#snprslist=fread('../genotypes/600K/Biogemma_LD_0.2_filtered_SNPs.txt',data.table=F,header=F)
us=fread('pgs/marker_effect_sizes.txt',data.table=F)
es=fread('pgs/Castelletti_DTA_marker_effect_sizes.txt',data.table=F)

us=us[us$marker %in% es$SNP.Name,]
m=intersect(es$SNP.Name, all_genos$marker)

all_genos=all_genos[all_genos$marker %in% m,]
rownames(all_genos)=seq(1,nrow(all_genos))

#keptsnps=fread('genotypes/600K/Biogemma_LD_0.2_filtered_SNPs.txt',data.table=F,header=F)
#genos=fread('genotypes/qtl2/Biogemma_DHgenos/rrBLUP_genotype.txt',data.table=F)
#all_genos=all_genos[all_genos$marker %in% snprslist$V1,]
#allres=c()
#gv_list=list()

#geno=all_genos[all_genos$chr==chr,]
#subu=us[us$chr==chr,]
X=all_genos[,sample.id]
  #colnames(X)=i
  #rownames(X)=geno$marker
rownames(X)=all_genos$marker
X=t(X)
X=as.matrix(X)
  #index=names(which(duplicated(X, MARGIN = 1)))
  #K=K[rownames(K)!=index,colnames(K)!=index]
  #y=data_blup[data_blup$ID!=index,c('y')]
  #X=[rownames(X)!=index,]
  #Zs <- scale(X, center = TRUE, scale = TRUE)
i=intersect(colnames(X),es$SNP.Name)
ms=colnames(X)
u=es[match(i,es$SNP.Name),]$effect_size
u=as.matrix(u)
  #Zs=Zs[,i]
#X=X[,i]
res=X %*% u
#allres=cbind(allres,res)

u2=us[match(i,us$marker),]$effect
u2=as.matrix(u2)
  #Zs=Zs[,i]
#X=X[,i]
res2=X %*% u2

#res=rowSums(allres)
res=as.data.frame(res,stringsAsFactors=F)
names(res)=c('dta')
res$ind=rownames(res)
res$rep=rep
fwrite(res,sprintf('pgs/Rep%.0f_Castelletti_DTA_FT_PGS.txt',rep),row.names=F,quote=F,sep='\t')

res2=as.data.frame(res2,stringsAsFactors=F)
names(res2)=c('dta')
res2$ind=rownames(res2)
res2$rep=rep
fwrite(res2,sprintf('pgs/Rep%.0f_Castelletti_markers_BG_DTA_FT_PGS.txt',rep),row.names=F,quote=F,sep='\t')
