#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('ggplot2')

# Read in Kinship Matrix
K=fread(sprintf('K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
# remove rows with missing phenotype data 
data=data[!is.na(data$y),]

m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

gwas_null=vector("list",length=16)
candidate_marker=""
min_pval=1
min_pval_hapgrp=0

# Run GridLMM
null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F,verbose=F)

h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup

Y=as.matrix(data_blup$y)
X_cov_null=null_model$lmod$X


for(h in seq(2,16)){
      X_list=readRDS(sprintf('../haplotype_probs/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))

      # Run GridLMM
      null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F,verbose=F)

      h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
      names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
      h2_start
      V_setup=null_model$setup

      Y=as.matrix(data_blup$y)
      X_cov=null_model$lmod$X
      X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,])
      
      X_list_null=NULL

      gwas=run_GridLMM_GWAS(Y,X_cov_null,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

      # If the variance of the haplotype probabilities is less than 1e-4, set the effect size to 0.

      #saveRDS(gwas,sprintf('models/Biogemma_chr%s_haplogrp%.0f_%s_x_ALL.rds',chr,h,pheno))
      sig=min(gwas$p_value_ML)
      if(sig < min_pval){
      	     min_pval=sig
	     min_pval_hapgrp=h
	     candidate_marker=gwas[which.min(gwas$p_value_ML),]$X_ID
      }
      i=h-1
      gwas_null[[i]]=gwas
}

threshtable=fread('threshold_table.txt',data.table=F)
cutoff=threshtable[threshtable$phenotype==pheno & threshtable$method=="haplotype_probs" & threshtable$environment=="ALL",]$threshold

print(candidate_marker)

#Create the candidate marker covariate matrix
X_list=readRDS(sprintf('../haplotype_probs/bg%s_filtered_haplogroup%.0f_probs.rds',chr,min_pval_hapgrp))
X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,])
null_columns=colnames(X_cov_null)
tmp=lapply(X_list_ordered,function(x) x[,candidate_marker])
tmp=as.data.frame(tmp,stringsAsFactors=F)
X_cov=cbind(X_cov_null,tmp)

names(X_cov)=c(null_columns,sprintf('hapgrp%.0f',seq(1,min_pval_hapgrp)))
rownames(X_cov)=rownames(X_cov_null)
X_cov=as.matrix(X_cov)

gwas_marker=vector("list",length=16)
for(h in seq(2,16)){
      X_list=readRDS(sprintf('../haplotype_probs/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
      X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,])

      if(h==min_pval_hapgrp){
	names=dimnames(X_list_ordered[[1]])[[2]]
	isit=names==candidate_marker

	location=which(isit==T,isit)
	X_minus=lapply(X_list_ordered,function(x) x[,-location])
	
      }else{
	X_minus=X_list_ordered
      }

      gwas_m=run_GridLMM_GWAS(Y,X_cov,X_minus[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method="ML",mc.cores=cores,verbose=F)
      i=h-1
      gwas_marker[[i]]=gwas_m
}

all_chrom=c()


for(h in seq(2,16)){
      i=h-1
      gwas_m=gwas_marker[[i]]
      gwas=gwas_null[[i]]

      if(h == min_pval_hapgrp){
      	   gwas_null_drop=gwas[gwas$X_ID!=candidate_marker,]
      }else{
	   gwas_null_drop=gwas
      }

      pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
      p=match(gwas_null_drop$X_ID,pmap$marker)
      phy_pos=pmap[p,]$pos
      chrom=data.frame(chrm=chr,hapgrp=h,ID=gwas_null_drop$X_ID,pos=phy_pos,null_pvalues=gwas_null_drop$p_value_ML,m_pvalues=gwas_m$p_value_ML,null_ML_logLik=gwas_null_drop$ML_logLik,m_ML_logLik=gwas_m$ML_logLik,stringsAsFactors=F)

      #min_pos=chrom[chrom$ID==candidate_marker,]$pos
      chrom$log10p=-log10(chrom$null_pvalues)
      chrom$sig = chrom$log10p >= cutoff
      rownames(chrom)=seq(1,dim(chrom)[1])

      chrom$lrt=(-log10(pchisq(chrom$m_ML_logLik-chrom$null_ML_logLik,1,lower.tail=F)))
      chrom=chrom[!is.na(chrom$null_pvalues),]
      chrom=chrom[!is.na(chrom$pos),]
      chrom=chrom[!is.na(chrom$m_pvalues),]

      all_chrom=rbind(all_chrom,chrom)
}

fwrite(all_chrom,sprintf('bound_images/%s_BLUP_haplotypeprobs_chr%s_LRT.txt',pheno,chr),row.names=F,quote=F)

png(sprintf('bound_images/%s_BLUP_haplotypeprobs_chr%s_LRT.png',pheno,chr),width=800,height=800)
print(ggplot(all_chrom,aes(x=pos,y=lrt)) + geom_point(aes(color=sig),alpha=0.5) + scale_color_manual(breaks=chrom$sig,values=c("FALSE"="black","TRUE"="red")) + geom_hline(yintercept=cutoff,color="green") + geom_hline(yintercept=(-log10(min_pval)),color="blue") + guides(color=F) + xlab("Position (Mb)") + ylab("-log10(p-value)") + ggtitle(sprintf("%s BLUP Haplotype Prob LRT Chr %s",pheno,chr)) + theme_classic())
dev.off()

png(sprintf('bound_images/%s_BLUP_haplotypeprobs_chr%s_LRT_by_hapgrp.png',pheno,chr),width=800,height=800)
print(ggplot(all_chrom,aes(x=pos,y=lrt)) + geom_point(aes(color=as.factor(hapgrp)),alpha=0.5) + geom_hline(yintercept=cutoff,color="green") + geom_hline(yintercept=(-log10(min_pval)),color="blue") + xlab("Position (Mb)") + ylab("-log10(p-value)") + ggtitle(sprintf("%s BLUP Haplotype Prob LRT Chr %s",pheno,chr)) + theme_classic())
dev.off()