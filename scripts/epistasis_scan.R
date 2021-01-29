#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('GridLMM')
library('data.table')
library('lme4')
library('reshape2')
library('tibble')
library('dplyr')
library('qtl2')
library('ggplot2')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite=c(T,F,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

plotdf_all=c()
tests=0
#for(c in 1:10){
#  print(c)
#  chr=as.character(c)

geno=fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
rownames(geno)=geno$ind
  #bg_clean=clean_genoprob(bg,column_threshold=0.01,value_threshold = 0.01,cores=1)

mite_start=135.947816
mite_end=135.946644
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

inds=geno$ind
  # Read in phenotype BLUP data for flowering time
phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
m1=lmer(male_flowering_d6~Loc.Year.Treat + (1|Genotype_code),phenotype)
data_blup = as.data.frame(ranef(m1)$Genotype_code)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]
geno = geno[rownames(geno) %in% data_blup$ID,]
inds=inds[inds %in% data_blup$ID]
mat=match(inds,data_blup$ID)
df=data.frame(ind=inds,y=data_blup[mat,]$y,stringsAsFactors=F)

mite=fread('GridLMM/mite_probabilities.txt',data.table=F)
mite=mite[mite$ID %in% data_blup$ID,]
rownames(mite)=seq(1,nrow(mite))
mite$clean_m=ifelse(mite$`AX-91102970`>=0.9,1,0)
mite_marker=names(mite)[2]

markers=names(geno)[-1]
  #markers=markers[-41456]
if(chr=="8"){
  geno=geno[,-(which(names(geno)==mite_marker))]
}
markers=names(geno)[-1]
tots=colSums(geno[,markers])
  #markers=markers[!(markers %in% names(which(tots<=15)))]
results=list()
count=0

for(m in markers){
  subgeno=geno[,c('ind',m)]
  subdf=df
  subdf$marker=geno[match(subdf$ind,geno$ind),m]
  subdf$mite=mite[match(subdf$ind,mite$ID),'clean_m']
  names(subdf)=c('ind','y','test_marker','mite')
  subdf$test=paste0(subdf$test_marker,'_',subdf$mite)
  if(length(unique(subdf$test)) == 4){
    model=lm(y ~ as.factor(test_marker)*as.factor(mite), data=subdf)
    count=count+1
    results[[count]]=list(pvalue=summary(model)$coeff[[16]],test_marker=m)
  }
  else{
    count=count+1
    results[[count]]=list(pvalue=NA,test_marker=m)
  }
}
pvalues=unlist(lapply(results,function(x) -log10(x$pvalue)))
markers=unlist(lapply(results,function(x) x$test_marker))
plotdf=data.frame(marker=markers,log10p=pvalues,stringsAsFactors=F)
plotdf$pos=pmap[match(plotdf$marker,pmap$marker),]$pos
plotdf=plotdf[!is.na(plotdf$log10p),]
plotdf$chrom==as.numeric(chr)
  #tests=tests+nrow(plotdf)
  #saveRDS(results,sprintf('GridLMM/MITE_only/chr%s_vgt1_epistasis_scan.rds',chr))
  #plotdf_all=rbind(plotdf_all,plotdf)



#cutoff=-log10(0.05/tests)

fwrite(plotdf,sprintf('GridLMM/MITE_only/chr%s_vgt1_epistasis_scan.txt',chr),quote=F,sep='\t',row.names=F)

#png('GridLMM/MITE_only/vgt1_epistasis_scan.png',width=1000,height=800)
#print(ggplot(plotdf,aes(x=pos/1e6,y=log10p)) + geom_point() + geom_hline(yintercept=cutoff) + facet_grid(.~chrom))
#dev.off()
