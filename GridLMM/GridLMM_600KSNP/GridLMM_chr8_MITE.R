#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
cores=as.numeric(args[[1]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('reshape2')
library('tibble')
library('dplyr')
library('ggplot2')

pheno="male_flowering_d6"
chr="8"
founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite=c(T,F,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

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

#map=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

X=fread(sprintf('../../genotypes/qtl2/Biogemma_DHgenos/bg%s_filtered_600K.csv',chr),data.table=F,stringsAsFactors=F)
#X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
#X=t(X)
#print(dim(X))
#X=as.data.frame(X)
#rownames(X)=geno$ind
#colnames(X)=colnames(geno)[2:dim(geno)[2]]
rownames(X)=X$ind
X=X[,2:dim(X)[2]]

X=X[rownames(X) %in% data_blup$ID,]
X=as.matrix(X)
data_blup=data_blup[data_blup$ID %in% rownames(X),]

dimr=dim(X)[1]
#dimc=dim(X)[2]
mono=c()
m=apply(X,MARGIN=2,FUN=function(n) length(unique(n)))
if(length(m[m==dimr]!=0)){
  mono=c(mono,which(m==dimr))
}
if(length(mono)>=1){
  X_filtered=X[,-mono]
}else{
  X_filtered=X
}

X=as.matrix(X_filtered)
marker='AX-91102970'
at_mite=fread('../mite_probabilities.txt',data.table=F)
at_mite=at_mite[at_mite$ID %in% data_blup$ID,]
rownames(at_mite)=at_mite$ID
at_mite=unlist(unname(at_mite[data_blup$ID,marker]))
data_blup$at_mite=at_mite
# Run GridLMM
gwas = GridLMM_GWAS(
                        formula = y~at_mite + (1|ID),
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
                        centerX = TRUE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'ML',
                        mc.cores = cores,
                        verbose = FALSE
)

saveRDS(gwas,sprintf('models/chr%s_%s_x_BLUP_MITE_covariate_600KSNP_ML.rds',chr,pheno))

mod=readRDS(sprintf('models/chr%s_%s_x_BLUP_MITE_covariate_600KSNP_ML.rds',chr,pheno))
p=match(mod$results$X_ID,pmap$marker)
phy_pos=pmap[p,]$pos
ID=mod$results$X_ID
gwas_results=data.frame(chr=chr,pos=phy_pos,pvalues=mod$results$p_value_ML,stringsAsFactors=F)
gwas_results$log10p = -log10(gwas_results$pvalues)

thresh_table=fread('../threshold_table.txt',data.table=F,stringsAsFactors=F)
rec=thresh_table$phenotype==pheno & thresh_table$environment=="ALL" & thresh_table$method=="600K_SNP"
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
 ggtitle(sprintf("DA  Chrom %s with MITE as Covariate Using 600K SNPs",chr,pheno)) +
  xlab("Distance (Mb)") + ylab("-log10(P-Value)") +
  geom_vline(xintercept=mite_start,color="black") +
   guides(color=F) + labs(caption=label))
dev.off()

og_mod=readRDS(sprintf('models/chr%s_%s_x_ALL_600KSNP_ML.rds',chr,pheno))
og_mod=og_mod$results
mod=mod$results
sig_markers=unique(c(og_mod[-log10(og_mod$p_value_ML)>=cutoff,]$X_ID,mod[-log10(mod$p_value_ML)>=cutoff,]$X_ID))

comp=data.frame(marker=sig_markers)
comp$null=-log10(og_mod[match(sig_markers,og_mod$X_ID),]$p_value_ML)
comp$covariate=-log10(mod[match(sig_markers,mod$X_ID),]$p_value_ML)
comp$null_sig=comp$null>=cutoff
comp$cov_sig=comp$covariate>=cutoff
comp$pos=pmap[match(sig_markers,pmap$marker),]$pos

fwrite(comp,'MITE_model_comparison_sig_SNPs.txt',row.names=F,quote=F,sep='\t')
