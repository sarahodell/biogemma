#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('ggplot2')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")

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
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y = data$y - mean(data$y)
rownames(data)=seq(1,dim(data)[1])

# Read in the haplotype group probabilities
# Filter genotypes that are not in the K matrix
X_list=readRDS(sprintf('../../geno_probs/bg%s_filtered_genotype_probs.rds',chr))

# Run GridLMM
null_model = GridLMM_ML(y~1 + (1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup

Y=as.matrix(data$y)
X_cov_null=null_model$lmod$X
#head(X_cov_null)

X_list_ordered=lapply(X_list,function(x) x[data$ID,])
X_list_null=NULL

gwas_null=run_GridLMM_GWAS(Y,X_cov_null,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method="ML",mc.cores=cores,verbose=F)


threshtable=fread('../threshold_table.txt',data.table=F)
cutoff=threshtable[threshtable$phenotype==pheno & threshtable$method=="founder_probs" & threshtable$environment==env,]$threshold

candidate_marker=gwas_null[which.min(gwas_null$p_value_ML),]$X_ID
min_pval=gwas_null[which.min(gwas_null$p_value_ML),]$p_value_ML

names=dimnames(X_list_ordered[[1]])[[2]]
isit=names==candidate_marker

location=which(isit==T,isit)
X_minus=lapply(X_list_ordered,function(x) x[,-location])

X_cov=X_cov_null
null_columns=colnames(X_cov_null)
tmp=lapply(X_list_ordered, function(x) x[data$ID,location])
tmp=as.data.frame(tmp,stringsAsFactors=F)
X_cov=cbind(X_cov,tmp)
names(X_cov)=c(null_columns,founders)
rownames(X_cov)=rownames(X_cov_null)

#head(X_cov)
#dim(X_cov)
#is(X_cov)
X_cov=as.matrix(X_cov)

gwas_marker=run_GridLMM_GWAS(Y,X_cov,X_minus[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method="ML",mc.cores=cores,verbose=F)

gwas_null_drop=gwas_null[gwas_null$X_ID!=candidate_marker,]

pmap=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
p=match(gwas_null$X_ID,pmap$marker)
phy_pos=pmap[p,]$pos
chrom=data.frame(chrm=chr,ID=gwas_null$X_ID,pos=phy_pos,pvalues=gwas_null$p_value_ML,stringsAsFactors = F)

min_pos=chrom[chrom$ID==candidate_marker,]$pos

chrom$log10p=-log10(chrom$pvalues)
chrom$sig = chrom$log10p >= cutoff
rownames(chrom)=seq(1,dim(chrom)[1])

chrom=chrom[chrom$ID!=candidate_marker,]
chrom$lrt=(-log10(pchisq(gwas_marker$ML_logLik-gwas_null_drop$ML_logLik,1,lower.tail=F)))
chrom=chrom[!is.na(chrom$pvalues),]
chrom=chrom[!is.na(chrom$pos),]

fwrite(chrom,sprintf('bound_images/%s_%s_founderprobs_chr%s_LRT.txt',pheno,env,chr),row.names=F,quote=F)

png(sprintf('bound_images/%s_%s_founderprobs_chr%s_LRT.png',pheno,env,chr),width=800,height=800)
print(ggplot(chrom,aes(x=pos,y=lrt)) + geom_point(aes(color=sig),alpha=0.5) + scale_color_manual(breaks=chrom$sig,values=c("FALSE"="black","TRUE"="red")) + geom_hline(yintercept=cutoff,color="green") + geom_hline(yintercept=(-log10(min_pval)),color="blue") + guides(color=F)+ xlab("Position (Mb)") + ylab("-\
log10(p-value)") + ggtitle(sprintf("%s %s Founder Prob LRT Chr %s",pheno,env,chr)) +theme_classic())
dev.off()



#saveRDS(gwas,sprintf('models/Biogemma_chr%s_%s_x_ALL_founderprobs.rds',chr,pheno))

# Convert all very high and very low probabilities to 1 and 0, respectively
#X_list_full = lapply(X_list_ordered,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.95,1,ifelse(x[,i]<=0.05,0,x[,i]))))
#for(i in 1:16){dimnames(X_list_full[[i]])[[2]]=dimnames(X_list_ordered[[i]])[[2]]}
#
#gwas_adjusted=gwas
#sums=lapply(X_list_full,function(x) colSums(x))
#for(i in 1:16){
#    s=sums[[i]]
#    t=dim(X_list_full[[i]])[1]-2
#    l=2
#    grab=which(s>t,s)
#    grab=c(grab,which(s<l,s))
#    grab=sort(grab)
#    beta=sprintf('beta.%.0f',seq(1,16))
#    gwas_adjusted[grab,beta]=0
#    print(grab)
#}

#saveRDS(gwas_adjusted,sprintf('models/Biogemma_chr%s_%s_x_ALL_founderprobs_adjusted.rds',chr,pheno))
