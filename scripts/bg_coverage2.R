#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

### Analyzing DH line genotype probabilities


library("qtl2")
library("dplyr")
library("data.table")
library("ggplot2")
library("reshape2")
library("tibble")


all_pr=readRDS(sprintf('bg%s_genoprobs_010319.rds',c))
all_pr_c=clean_genoprob(all_pr)
pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

hex_colors=c("#f42896","#84ef7c","#234489","#8ed1d6","#702349","#f2875b","#f1ff00","#bbbbbb","#ff2a3a","#56cc59","#663dd3","#478959","#47145\
b",'#0f0e0e',"#ad147a", "#afb735","#ff5a00","#fc1919")
founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra",
"FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

coverage=c()

for(f in 1:16){
    df=all_pr_c[[1]][,f,]
    tm =  t(df)
    mdf=as.data.frame(tm,row.names=rownames(tm))
    mdf=rownames_to_column(mdf,"marker")
    mdf = merge(mdf,pmap,by.x="marker",by.y="marker")
    mdf = mdf[,c(2:345,347)]
    mlong=melt(mdf,id="pos")
    # take only sites where there is greater than 0.95 probability of coming from that founder
    mlong=mlong[mlong$value>0.95,]
    presence=mlong %>% group_by(pos) %>% summarize(total=length(value))
    presence=as.data.frame(presence)
    max_prob=data.frame(pos=presence$pos,total=presence$total,founder=rep(founders[f],dim(presence)[1]),chr=rep(c,dim(presence)[1]))
    coverage=rbind(coverage,max_prob)
}
cov_counts=coverage %>% group_by(pos,founder) %>% summarize(tot_cov=round(total,2))

fwrite(cov_counts,sprintf('founder_coverage_counts_chr%s.txt',c),sep='\t',row.names=F,quote=F)

png(sprintf('founder_counts_chr%s.png',c),width=960,height=960)
print(ggplot(cov_counts,aes(x=pos/1e6,y=tot_cov)) + geom_point(aes(color=founder)) + facet_grid(founder~.) + ggtitle(sprintf("Founder Coverage Chromosome %s",c)) + xlab("Position (Mb)") + ylab("Counts") + geom_hline(yintercept=21.5,color="black") + theme_classic())
dev.off()
