#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])
pheno=as.character(args[[2]])
env=as.character(args[[3]])

library('data.table')
library('ggplot2')
library('dplyr')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
baselist=c(8,7,7,8,6,7,7,8,7,7)

ibd=fread(sprintf('../../ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',c),data.table=F)
pmap=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

all_effects=c()

for(h in baselist[as.numeric(c)]:16){
      model=readRDS(sprintf('models/Biogemma_chr%s_haplogrp%.0f_%s_x_%s.rds',c,h,pheno,env))
      model_merge=merge(model,pmap,by.x="X_ID",by.y="marker")
      model_merge=model_merge[order(model_merge$pos),]
      row.names(model_merge)=seq(1,dim(model_merge)[1])
      for(f in founders){
      	    effects=c()
      	    sub_ibd=ibd[ibd$n_haps==h,c('chrom','start','end',f)]
      	    sub_ibd$betas=sprintf('beta.%.0f',sub_ibd[,f])
      	    for(i in 1:dim(model_merge)[1]){
      	    	  pos=model_merge[i,'pos']
      	    	  within=sub_ibd[(pos>=sub_ibd$start) & (pos<sub_ibd$end),]$betas
      	    	  eff=model_merge[i,within]
      	    	  effects=c(effects,eff)
      	    }
	    tmp=data.frame(chr=c,hapgrp=h,pos=model_merge$pos,marker=model_merge$X_ID,effect_size=effects,founder=f)
	    all_effects=rbind(all_effects,tmp)
      }
}

all_effects=as.data.frame(all_effects)
names=c('chr','hapgrp','pos','marker','effect_size','founder')
all_effects$hapgrp=as.factor(all_effects$hapgrp)

png(sprintf('images/chr%s_%s_x_%s_effect_sizes_by_founder_by_hapgrp.png',c,pheno,env),width=960,height=960)
print(ggplot(all_effects,aes(x=pos/1e6,y=effect_size)) + geom_point(aes(color=hapgrp)) + facet_wrap(~founder) + theme_classic() + theme(strip.text.x = element_text(size = 16)))
dev.off()

bounds=fread('../Biogemma_QTL.csv',data.table=F)
mite_start=135.947816
mite_end=135.946644

left_bound=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$left_bound_bp
right_bound=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$alt_right_bound_bp

if(left_bound/1e6 > mite_start){
  left_bound=135e6
}
sub8=all_effects[all_effects$pos>left_bound & all_effects$pos < right_bound,]

png(sprintf('images/chr%s_%s_x_%s_effect_sizes_by_founder_by_hapgrp_zoom_vgt1.png',c,pheno,env),width=960,height=960)
print(ggplot(sub8,aes(x=pos/1e6,y=effect_size)) + geom_point(aes(color=hapgrp)) + facet_wrap(~founder) + geom_vline(xintercept=135.9) + geom_hline(yintercept=0,color="red") + ggtitle("Chr 8:107-173 Founder Effect Sizes") + theme_classic() + theme(strip.text.x = element_text(size = 16)))
dev.off()

#png(sprintf('chr%s_log10_effect_sizes_by_founder_by_hapgrp.png',c),width=960,height=960)
#print(ggplot(all_effects,aes(x=pos/1e6,y=log10(abs(effect_size)))) + geom_point(aes(color=hapgrp)) + facet_wrap(~founder))
#dev.off()

png(sprintf('images/chr%s_%s_x_%s_effect_sizes_by_founder.png',c,pheno,env),width=800,height=530)
print(ggplot(all_effects,aes(x=pos/1e6,y=effect_size)) + geom_jitter(aes(color=founder)) + geom_vline(xintercept=135) + ggtitle(sprintf("%s %s Effect Sizes",pheno,env)) + labs(subtitle="Chr. 8, Haplotype Probabilities") + theme(strip.text.x = element_text(size = 16)) + xlab("Position (Mb)") + ylab("Effect Size (ggd)"))
dev.off()
