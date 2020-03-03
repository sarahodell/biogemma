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


#ibd=fread(sprintf('../../ibd_segments/bg%s_ibd_blocks.txt',c),data.table=F)
pmap=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

all_effects=c()

model=readRDS(sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs_adjusted.rds',c,pheno,env))
model_merge=merge(model,pmap,by.x="X_ID",by.y="marker")
model_merge=model_merge[order(model_merge$pos),]
row.names(model_merge)=seq(1,dim(model_merge)[1])
betas=sprintf('beta.%.0f',seq(1,16))
sub = model_merge[,c(betas)]
names(sub)=founders
sub$pos=model_merge$pos
mlong=melt(sub,'pos')

png(sprintf('chr%s_%s_x_%s_founderprobs_effect_sizes_by_by_founder.png',c,pheno,env),width=960,height=960)
print(ggplot(mlong,aes(x=pos,y=value)) + geom_point(aes(color=variable)) + facet_wrap(~variable) + geom_hline(yintercept=0,color="red") + guides(color=F) + ggtitle(sprintf("Founder Probability %s Effect Sizes in %s, Chr.%s",pheno,env,c)) + ylab("Effect Size") + xlab("Position (Mb)"))
dev.off()