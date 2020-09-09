#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])
pheno=as.character(args[[2]])
env=as.character(args[[3]])

library('data.table')
library('ggplot2')
library('dplyr')
library('reshape2')
#library('emmeans')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")


#ibd=fread(sprintf('../../ibd_segments/bg%s_ibd_blocks.txt',c),data.table=F)
pmap=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

all_effects=c()

model=readRDS(sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs.rds',c,pheno,env))
model_merge=merge(model,pmap,by.x="X_ID",by.y="marker")
model_merge=model_merge[order(model_merge$pos),]
row.names(model_merge)=seq(1,dim(model_merge)[1])
betas=sprintf('beta.%.0f',seq(1,16))
sub = model_merge[,c(betas)]
names(sub)=founders
sub$pos=model_merge$pos
mlong=reshape2::melt(sub,'pos')

bounds=fread('../Biogemma_QTL.csv',data.table=F)
mite_start=135.947816
mite_end=135.946644

left_bound=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$left_bound_bp
right_bound=bounds[bounds$Chromosome==c & bounds$Phenotype==pheno & bounds$Method=="Founder_probs" & bounds$Environment==env,]$alt_right_bound_bp

if(left_bound/1e6 > mite_start){
  left_bound=135e6
}

mlong_sub=mlong[mlong$pos>left_bound & mlong$pos<right_bound,]

png(sprintf('images/chr%s_%s_x_%s_founderprobs_effect_sizes_by_by_founder.png',c,pheno,env),width=960,height=960)
print(ggplot(mlong_sub,aes(x=pos/1e6,y=value)) + geom_point(aes(color=variable)) + facet_wrap(~variable) + geom_hline(yintercept=0,color="red") +
geom_vline(xintercept=mite_start,color="black") +
 guides(color=F) + ggtitle(sprintf("Founder Probability %s Effect Sizes in %s, Chr.%s:%.2f-%.2f Mb",pheno,env,c,left_bound/1e6,right_bound/1e6)) + ylab("Effect Size") + xlab("Position (Mb)"))
dev.off()
