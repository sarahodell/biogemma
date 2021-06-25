#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.numeric(args[[1]])

#library('LDheatmap')
library('data.table')
library('tidyverse')
library('ggplot2')
library('cowplot')


#awk -F';' '{print > "Biogemma_DHLines_chr"$1"_rsquared.ld.txt"}' Biogemma_DHLines_rsquared.ld
#ld=fread('/scratch/sodell/Biogemma_DHLines_rsquared.ld',data.table=F)
ld=fread(sprintf('stats/ld_decay/%.0f_rsquared.ld',chr),data.table=F)
#map=fread('genotypes/plink_files/600K/Biogemma_Founders_600K.map',data.table=F)
names(ld)=c('CHR_A','BP_A','SNP_A','MAF_A','CHR_B','BP_B','SNP_B','MAF_B','R2')

#ld=ld[ld$CHR_A==chr & ld$CHR_B == chr,]
#fwrite(ld,sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_rsquared.ld',chr),row.names=F,quote=F,sep='\t')
f_freq=fread('stats/ld_decay/Biogemma_Founders_allele_freq.frqx',data.table=F)

ld$F_MAF_A=f_freq[match(ld$SNP_A,f_freq$SNP),]$MAF
ld$F_MAF_B=f_freq[match(ld$SNP_B,f_freq$SNP),]$MAF
ld$dist=ld$BP_B - ld$BP_A
remove(f_freq)
#100kb bins
bins=seq(1,max(ld$dist)+1000,100000)
ld$group=cut(ld$dist,breaks=bins,labels=seq(1,length(bins)-1))

decay=ld %>% group_by(group) %>% summarize(r2=mean(R2))
rm(ld)
decay$chr=chr
decay$group2=as.numeric(decay$group)
decay=decay[order(decay$group2),]
#rownames(decay)=seq(1,nrow(decay))
fwrite(decay,sprintf('stats/ld_decay/chr%.0f_mean_R2.txt',chr),row.names=F,quote=F,sep='\t')

png(sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_ld_decay.png',chr),width=1000,height=1000)
print(ggplot(decay,aes(x=group2,y=r2,group=group2)) +  geom_point() + geom_line() + xlab("Distance (10kb)") + ylab("R-squared"))
dev.off()


plot_list=list()
for(c in 1:10){
  decay=fread(sprintf('stats/ld_decay/chr%.0f_mean_R2.txt',c),data.table=F)
  plot_list[[c]]=ggplot(decay,aes(x=group2,y=r2,group=group2)) +
    geom_point() + geom_line() + xlab("Distance (100kb)") +
    scale_y_continuous(limits=c(0,0.4),breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)) +
     ylab("R-squared")
}

pall=plot_grid(plotlist=plot_list,ncol=3,labels=c('1','2','3','4','5','6','7','8','9','10'))
png('stats/ld_decay/all_chroms_ld_decay.png',width=1000,height=1000)
print(pall)
dev.off()

#decay$group2=as.numeric(decay$group)
#mafs=unique(ld$MAF_A)
#models=list()
#counts=1
#for(m in mafs){
#  mod=lm(R2 ~ dist^2,ld[ld$MAF_A==m,])
#  models[[counts]]=mod
#}
#names(models)=mafs

#png(sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_ld_decay.png',chr),width=1000,height=1000)
#print(ggplot(ld,aes(x=dist,y=R2)) + geom_point(aes(color=F_MAF_B)) + geom_smooth(se=F,color="black") + facet_wrap(~F_MAF_A))
#dev.off()


#png(sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_ld_decay_MAF.png',chr),width=1000,height=1000)
#print(ggplot(ld,aes(x=dist,y=R2)) + geom_point(aes(color=F_MAF_B)) + geom_smooth(se=F,color="black") + facet_grid(F_MAF_A~F_MAF_B))
#dev.off()

#ld_matrix = ld[,c('SNP_A','SNP_B','R2')] %>%
#  group_by(`SNP_A`) %>%
#  mutate(id = row_number()) %>%
#pivot_wider(names_from='SNP_A',values_from='R2')
#ld_matrix=as.data.frame(ld_matrix,stringsAsFactors=F)
#rownames(ld_matrix)=ld_matrix$SNP_B
#map=map[map$V2 %in% rownames(ld_matrix),]
#ld_matrix=ld_matrix[,-1]
#ld_matrix=as.matrix(ld_matrix)


#png('stats/ld_decay/Biogemma_Founders_LD_heatmap.png')
#LDheatmap(ld_matrix,genetic.distances=map$V3,distances=map$V4)
#dev.off()
