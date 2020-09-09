#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])
pheno=as.character(args[[2]])

library('data.table')
library('ggplot2')
library('dplyr')
library('reshape2')
#library('emmeans')
environments=c("ALL","BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD","STPAUL_2017_WD","SZEGED_2017_OPT")
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")


#ibd=fread(sprintf('../../ibd_segments/bg%s_ibd_blocks.txt',c),data.table=F)
pmap=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)
bounds=fread('../Biogemma_QTL.csv',data.table=F)
mlong_all=c()
mlong_point=c()
mite_start=135.947816
mite_end=135.946644

for(env in environments){
  model=readRDS(sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs.rds',c,pheno,env))
  model_merge=merge(model,pmap,by.x="X_ID",by.y="marker")
  model_merge=model_merge[order(model_merge$pos),]
  row.names(model_merge)=seq(1,dim(model_merge)[1])
  betas=sprintf('beta.%.0f',seq(1,16))
  sub = model_merge[,c(betas)]
  names(sub)=founders
  sub$pos=model_merge$pos
  mlong=melt(sub,'pos')
  mlong$env=env
  left_bound=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$left_bound_bp
  right_bound=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$alt_right_bound_bp
  highest_snp=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$highest_SNP
  if(left_bound/1e6 > mite_start){
    left_bound=135e6
  }
  highest_pos=pmap[pmap$marker==highest_snp,]$pos
  mlong_sub=mlong[mlong$pos>left_bound & mlong$pos<right_bound,]
  mlong_snp=mlong[mlong$pos==highest_pos,]
  mlong_all=rbind(mlong_all,mlong_sub)
  mlong_point=rbind(mlong_point,mlong_snp)
}
mlong_all=as.data.frame(mlong_all,stringsAsFactors=F)
names(mlong_all)=c("pos","variable","value","env")


png(sprintf('images/chr%s_%s__founderprobs_avg_effect_sizes_all_founders.png',c,pheno,env),width=960,height=960)
print(ggplot(mlong_all,aes(x=factor(env),y=value)) + geom_boxplot(aes(fill=env)) + facet_wrap(~variable) +
ggtitle(sprintf("Founder Probability %s  Average vgt1 Effect Sizes by Environment (within vgt1 QTL bound for chromosome %s)",pheno,c)) + ylab("Effect Size") + geom_hline(yintercept=0,color='red') +
 xlab("Environment") +
 theme_classic() +
 theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major.y = element_line(color = "grey80"),panel.grid.minor.y = element_line(color = "grey80")))
dev.off()

png(sprintf('images/chr%s_%s__founderprobs_highest_SNP_effect_sizes_all_founders.png',c,pheno,env),width=960,height=960)
print(ggplot(mlong_point,aes(x=factor(env),y=value)) + geom_bar(aes(fill=env),stat="identity") + facet_wrap(~variable) +
ggtitle(sprintf("Founder Probability %s vgt1 Effect Sizes by Environment at highest SNP, chromosome %s)",pheno,c)) + ylab("Effect Size") + geom_hline(yintercept=0,color='red') +
 xlab("Environment") +
 theme_classic() +
 theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major.y = element_line(color = "grey80"),panel.grid.minor.y = element_line(color = "grey80")))
dev.off()
