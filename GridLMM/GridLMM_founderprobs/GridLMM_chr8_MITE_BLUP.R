#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
cores=as.numeric(args[[1]])

#date=format(Sys.time(),'%m%d%y')

pheno="male_flowering_d6"
chr="8"

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('ggplot2')
library('reshape2')
library('tibble')
library('dplyr')

# Read in Kinship Matrix
K=fread(sprintf('../K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


pmap=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
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

# Read in the haplotype group probabilities
# Filter genotypes that are not in the K matrix
X_list=readRDS(sprintf('../../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

new_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
# Ordered for new_founders

#Make B73 the first in the list so that it is the one that is dropped
names(X_list)=founders
X_list=X_list[new_founders]

marker='AX-91102970'
at_mite=fread('../mite_probabilities.txt',data.table=F)
at_mite=at_mite[at_mite$ID %in% data_blup$ID,]
rownames(at_mite)=at_mite$ID
at_mite=unlist(unname(at_mite[data_blup$ID,marker]))
# Run GridLMM
null_model = GridLMM_ML(y~at_mite + (1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F,verbose=F)

h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup

Y=as.matrix(data_blup$y)
X_cov=null_model$lmod$X
X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,])
X_list_null=NULL

gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

saveRDS(gwas,sprintf('models/Biogemma_chr%s_%s_BLUP_MITE_covariate_founderprobs.rds',chr,pheno))

# Convert all very high and very low probabilities to 1 and 0, respectively
X_list_full = lapply(X_list_ordered,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.95,1,ifelse(x[,i]<=0.05,0,x[,i]))))
for(i in 1:16){dimnames(X_list_full[[i]])[[2]]=dimnames(X_list_ordered[[i]])[[2]]}

gwas_adjusted=gwas
sums=lapply(X_list_full,function(x) colSums(x))
for(i in 1:16){
    s=sums[[i]]
    t=dim(X_list_full[[i]])[1]-2
    l=2
    grab=which(s>t,s)
    grab=c(grab,which(s<l,s))
    grab=sort(grab)
    beta=sprintf('beta.%.0f',seq(1,16))
    gwas_adjusted[grab,beta]=0
    gwas_adjusted[grab,'p_value_ML']=0.99
    print(grab)
}

saveRDS(gwas_adjusted,sprintf('models/Biogemma_chr%s_%s_x_BLUP_MITE_covariate_founderprobs_adjusted.rds',chr,pheno))

#pmap=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_cs.csv',chr),data.table=F,stringsAsFactors=F)
mod=readRDS(sprintf('models/Biogemma_chr%s_%s_x_BLUP_MITE_covariate_founderprobs_adjusted.rds',chr,pheno))
p=match(mod$X_ID,pmap$marker)
phy_pos=pmap[p,]$pos
ID=mod$X_ID
gwas_results=data.frame(chr=chr,pos=phy_pos,pvalues=mod$p_value_ML,stringsAsFactors=F)
gwas_results$log10p = -log10(gwas_results$pvalues)

thresh_table=fread('../threshold_table.txt',data.table=F,stringsAsFactors=F)
rec=thresh_table$phenotype==pheno & thresh_table$environment=="ALL" & thresh_table$method=="founder_probs"
cutoff=thresh_table[rec,]$threshold
print(cutoff)

gwas_results$sig=gwas_results$log10p >= cutoff
label<-sprintf("5%% Permutation Significance Threshold = %.3f",cutoff)
mite_start=135.947816
mite_end=135.946644

png(sprintf('images/%s_x_BLUP_chr8_MITE_manhattan.png',pheno),width=960,height=680)
theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(size=26),axis.title=element_text(size=14,face="bold"))
theme_update(panel.background=element_blank())
print(ggplot(gwas_results,aes(x=pos/1e6,y=log10p)) +
geom_point(aes(color=sig)) +
 scale_color_manual(breaks=gwas$sig,values=c("FALSE"="black","TRUE"="red")) +
 ggtitle(sprintf("DA  Chrom %s with MITE as Covariate Using Founder Probabiliites",chr,pheno)) +
  xlab("Distance (Mb)") + ylab("-log10(P-Value)") +
  geom_vline(xintercept=mite_start,color="black") +
   guides(color=F) + labs(caption=label))
dev.off()

og_mod=readRDS(sprintf('models/Biogemma_chr%s_%s_x_ALL_founderprobs_adjusted.rds',chr,pheno))
sig_markers=unique(c(og_mod[-log10(og_mod$p_value_ML)>cutoff,]$X_ID,mod[-log10(mod$p_value_ML)>cutoff,]$X_ID))

comp=data.frame(marker=sig_markers)
comp$null=-log10(og_mod[match(sig_markers,og_mod$X_ID),]$p_value_ML)
comp$covariate=-log10(mod[match(sig_markers,mod$X_ID),]$p_value_ML)
comp$null_sig=comp$null>=cutoff
comp$cov_sig=comp$covariate>=cutoff
comp$pos=pmap[match(sig_markers,pmap$marker),]$pos

fwrite(comp,'MITE_model_comparison_sig_SNPs.txt',row.names=F,quote=F,sep='\t')
