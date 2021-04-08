#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.numeric(args[[1]])

#library('LDheatmap')
library('data.table')
library('tidyverse')
library('ggplot2')

#ld=fread('/scratch/sodell/Biogemma_DHLines_rsquared.ld',data.table=F)
ld=fread('stats/ld_decay/Biogemma_DHLines_rsquared.ld',data.table=F)
#map=fread('genotypes/plink_files/600K/Biogemma_Founders_600K.map',data.table=F)
ld=ld[ld$CHR_A==chr & ld$CHR_B == chr,]
fwrite(ld,sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_rsquared.ld',chr),row.names=F,quote=F,sep='\t')
f_freq=fread('stats/ld_decay/Biogemma_Founders_allele_freq.frqx',data.table=F)

ld$F_MAF_A=f_freq[match(ld$SNP_A,f_freq$SNP),]$MAF
ld$F_MAF_B=f_freq[match(ld$SNP_B,f_freq$SNP),]$MAF
ld$dist=ld$BP_B - ld$BP_A
remove(f_freq)

#mafs=unique(ld$MAF_A)
#models=list()
#counts=1
#for(m in mafs){
#  mod=lm(R2 ~ dist^2,ld[ld$MAF_A==m,])
#  models[[counts]]=mod
#}
#names(models)=mafs

png(sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_ld_decay.png',chr),width=1000,height=1000)
print(ggplot(ld,aes(x=dist,y=R2)) + geom_point(aes(color=F_MAF_B)) + geom_smooth(se=F,color="black") + facet_wrap(~F_MAF_A))
dev.off()


png(sprintf('stats/ld_decay/Biogemma_DHLines_chr%.0f_ld_decay_MAF.png',chr),width=1000,height=1000)
print(ggplot(ld,aes(x=dist,y=R2)) + geom_point(aes(color=F_MAF_B)) + geom_smooth(se=F,color="black") + facet_grid(F_MAF_A~F_MAF_B))
dev.off()

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
